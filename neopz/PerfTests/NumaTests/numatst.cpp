/**
 * @file
 * @brief Performance tests on NUMA architecture
 * @author Edson Borin
 * @since 2012
 */

#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include <cstdlib>

#include "pzbfilestream.h" // TPZBFileStream, TPZFileStream
#include "pzmd5stream.h"

#include "pzlog.h"

#include <fstream>
#include <string>

#ifdef LOG4CXX
static LoggerPtr loggerconverge(Logger::getLogger("pz.converge"));
static LoggerPtr logger(Logger::getLogger("main"));
#endif

#include "pzskylmat.h"

//#include "timing_analysis.h"
#include "arglib.h"
#include "run_stats_table.h"

#ifdef HAS_GETRUSAGE
#include <sys/resource.h> // getrusage
#endif

#ifdef USING_TBB
#include "tbb/task_scheduler_init.h"
using namespace tbb;
// If you have issues with: dyld: Library not loaded: libtbb.dylib
// try setting the LD path. Ex:
//   export DYLD_FALLBACK_LIBRARY_PATH=/Users/borin/Desktop/neopz/tbb40_297oss/lib/
#endif

void help(const char* prg)
{
    cout << "Compute the Decompose_LDLt method for the matrix" << endl;
    cout << endl;
    cout << "Usage: " << prg << "-if file [-v verbose_level] [-b] "
    << "[-tot_rdt rdt_file] [-op matrix_operation] [-h]" << endl << endl;
    cout << "matrix_operation:" << endl;
    cout << " 0: Decompose_LDLt()" << endl;
    cout << " 1: Decompose_LDLt2() -- deprecated (not working)" << endl;
    cout << " 2: Decompose_Cholesky()" << endl;
    clarg::arguments_descriptions(cout, "  ", "\n");
}

clarg::argString ifn("-ifn", "input matrix file name (use -bi to read from binary files)", "matrix.txt");
clarg::argInt affinity("-af", "affinity mode (0=no affinity, 1=heuristi 1)", 0);
clarg::argInt verb_level("-v", "verbosity level", 0);
int verbose = 0;
/* Verbose macro. */
#define VERBOSE(level,...) if (level <= verbose) cout << __VA_ARGS__

clarg::argInt mop("-op", "Matrix operation", 1);
clarg::argBool br("-br", "binary reference. Reference decomposed matrix file format == binary.", false);
clarg::argBool bi("-bi", "binary input. Input file format == binary.", false);
clarg::argBool bd("-bd", "binary dump. Dump file format == binary.", false);
clarg::argBool h("-h", "help message", false);
clarg::argBool copy_matrix_inside_thread("-cot", "copy on thread - copy matrix inside thread.", false);
clarg::argInt mstats("-mstats", "Matrix statistics vebosity level.", 0);
clarg::argInt maxcol("-maxcol", "Limit computation to max column (Use Resize(maxcol)).", 0);
clarg::argString gen_dm_sig("-gen_dm_md5", "generates MD5 signature for decomposed matrix into file.", "decomposed_matrix.md5");
clarg::argString chk_dm_sig("-chk_dm_md5", "compute MD5 signature for decomposed matrix and check against MD5 at file.", "decomposed_matrix.md5");
clarg::argString chk_dm_error("-chk_dm_error", "check the decomposed matrix error against a reference matrix. (use -br to read from binary files)", "ref_decomposed_matrix.txt");
clarg::argDouble error_tol("-error_tol", "error tolerance.", 1.e-12);
clarg::argString dump_dm("-dump_dm", "dump decomposed matrix. (use -bd for binary format)", "dump_matrix.txt");
clarg::argInt cholesky_blk("-chol_blk", "Cholesky blocking factor", 256);

/* Run statistics. */
RunStatsTable total_rst("-tot_rdt",
                        "Whole program (total) statistics raw data table");

clarg::argInt nmats("-nmats", "Number of matrizes to decompose simultaneously.", 1);

class FileStreamWrapper
{
public:
    FileStreamWrapper(bool b) : binary(b)
    {}
    ~FileStreamWrapper() {}
    
    void OpenWrite(const std::string& fn)
    {
        if (binary)
            bfs.OpenWrite(fn);
        else
            fs.OpenWrite(fn);
    }
    
    void OpenRead(const std::string& fn)
    {
        if (binary)
            bfs.OpenRead(fn);
        else
            fs.OpenRead(fn);
    }
    
    operator TPZStream&()
    {
        if (binary)
            return bfs;
        else
            return fs;
    }
    
protected:
    
    bool binary;
    TPZFileStream  fs;
    TPZBFileStream bfs;
};

#include <sched.h>     //sched_getcpu

std::vector< TPZSkylMatrix<REAL>* > matrices;

#ifdef USING_LIBNUMA
cpu_set_t dies_mask_array[8];
cpu_set_t mask_core0;
cpu_set_t mask_L20;
cpu_set_t mask_die0;
cpu_set_t mask_proc0;
cpu_set_t mask_oddcores;
cpu_set_t mask_evencores;


void print_mask(cpu_set_t* mask)
{
    for (int i=0; i<64; i++) {
        cout << CPU_ISSET(i, mask)?"1":"0";
    }
}

#endif


// 4 processors
// 2 dies / processor
// 4 L2 caches per die
// 2 cores per L2 cache

// FPU is shared among the 2 cores under the same L2 cache

// Best Assign policy
// # of threads : policy
//            1 : any core
//      2, 3, 4 : different processors
//      5 - 16  : different dies.
void setup_masks()
{
#define SET_RANGE(mskp,start,end) CPU_ZERO(mskp); \
for (int i=start; i<end; i++) CPU_SET(i,mskp)
    
#ifdef USING_LIBNUMA
    SET_RANGE(&dies_mask_array[0],0,8);
    SET_RANGE(&dies_mask_array[1],8,16);
    SET_RANGE(&dies_mask_array[2],16,24);
    SET_RANGE(&dies_mask_array[3],24,32);
    SET_RANGE(&dies_mask_array[4],32,40);
    SET_RANGE(&dies_mask_array[5],40,48);
    SET_RANGE(&dies_mask_array[6],48,56);
    SET_RANGE(&dies_mask_array[7],56,64);
    
    SET_RANGE(&mask_proc0,0,16);
    SET_RANGE(&mask_die0,0,8);
    SET_RANGE(&mask_L20,0,2);
    SET_RANGE(&mask_core0,0,1);
    
    CPU_ZERO(&mask_oddcores);
    CPU_ZERO(&mask_evencores);
    for(int i=0; i<64; i+=2) {
        CPU_SET(i,&mask_evencores);
        CPU_SET(i+1,&mask_oddcores);
    }
    
    if (verbose >= 1) {
        cout << "mask core 0     : "; print_mask(&mask_core0); cout << endl;
        cout << "mask core L2 0  : "; print_mask(&mask_L20); cout << endl;
        cout << "mask core die 0 : "; print_mask(&mask_die0); cout << endl;
        cout << "mask core proc 0: "; print_mask(&mask_proc0); cout << endl;
        cout << "mask evencores  : "; print_mask(&mask_evencores); cout << endl;
        cout << "mask oddcores   : "; print_mask(&mask_oddcores); cout << endl;
    }
#endif
}


//  CPU_SET(cpus[idx],&mask);
//int cpus[] = {0, 16, 32, 48, 8, 24, 40, 54};
void set_affinity(int af, int tidx)
{
#ifdef USING_LIBNUMA
    cpu_set_t* msk = NULL;
    switch (af) {
            
        case 1: {
            msk = dies_mask_array + (tidx%8);
            break;
        }
        case 2: {
            msk = &mask_proc0;
            break;
        }
        case 3: {
            msk = &mask_die0;
            break;
        }
        case 4: {
            msk = &mask_L20;
            break;
        }
        case 5: {
            msk = &mask_core0;
            break;
        }
        case 6: {
            msk = &mask_evencores;
            break;
        }
        case 7: {
            msk = &mask_oddcores;
            break;
        }
        case 0: // Do not set affinity
            return;
        default: {
            VERBOSE(2, "Warning: -af " << af
                    << " has not been defined. Not setting affinity");
            return;
        }
    }
    
    if (verbose >= 2) {
        cout << "Thread " << tidx << " affinity mask = ";
        print_mask(msk);
        cout << endl;
    }
    
    sched_setaffinity(0, sizeof(cpu_set_t), msk);
#endif
}

TPZSkylMatrix<REAL> matrix;

void init_decompose(int idx)
{
    //  cpu_set_t mask;
    //  CPU_ZERO(&mask);
    //  CPU_SET(cpus[idx],&mask);
    if (affinity.get_value() > 0) {
        set_affinity(affinity.get_value(), idx);
    }
    matrices[idx] = new TPZSkylMatrix<REAL>(matrix);
}

void compute_decompose(int idx)
{
    
    TPZSkylMatrix<REAL>* matrix = matrices[idx];
    
#define CASE_OP(opid,method)				\
case opid:						\
matrix->method;					\
break
    
    switch (mop.get_value()) {
            CASE_OP(0,Decompose_LDLt());
        case 1:
            std::cerr << "ERROR: deprecated operation -- decompose LDLt2 is no longer implemented." << std::endl;
            break;
            CASE_OP(2,Decompose_Cholesky());
            CASE_OP(3,Decompose_Cholesky_blk(cholesky_blk.get_value()));
        default:
            std::cerr << "ERROR: Invalid matrix operation type." << std::endl;
    }
}

#include <pthread.h>

class thread_timer_t
{
public:
    thread_timer_t() {}
    void start()
    {start_time = getms();}
    void stop()
    {stop_time = getms();}
    uint64_t get_start() {return start_time; }
    uint64_t get_stop() {return stop_time; }
    uint64_t get_elapsed() {return stop_time-start_time; }
private:
    uint64_t getms()
    {
        timeval t;
        gettimeofday(&t,NULL);
        return (t.tv_sec*1000) + (t.tv_usec/1000);
    }
    uint64_t start_time;
    uint64_t stop_time;
};

int nthreads_initialized;
int nthreads;
bool wait_for_all_init;
std::vector<thread_timer_t> thread_timer;
pthread_cond_t  cond=PTHREAD_COND_INITIALIZER;
pthread_cond_t  main_cond=PTHREAD_COND_INITIALIZER;
pthread_mutex_t glob_mutex=PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t main_mutex=PTHREAD_MUTEX_INITIALIZER;
bool run_parallel;

class synchronized_threads_t
{
public:
    synchronized_threads_t()
    {
        run_parallel=false;
        init_routine=NULL;
        parallel_routine=NULL;
        nthreads=0;
        nthreads_initialized = 0;
        wait_for_all_init = true;
    }
    
    void execute_n_threads(unsigned n,
                           void (*init_routine)(int),
                           void (*parallel_routine)(int));
    
    struct thread_arg_t
    {
        thread_arg_t(int t,void (*ir)(int), void (*pr)(int),
                     pthread_mutex_t* mt, pthread_cond_t* cd,
                     pthread_cond_t* mcd) :
        tid(t), init_routine(ir), parallel_routine(pr),
        glob_mutex(mt), cond(cd), main_cond(mcd)
        {}
        int tid;
        void (*init_routine)(int);
        void (*parallel_routine)(int);
        pthread_mutex_t* glob_mutex;
        pthread_cond_t* cond;
        pthread_cond_t* main_cond;
    };
    
private:
    
    /** Initialization routine. Called before the cond mutex. */
    void (*init_routine)(int);
    /** Parallel routine. Called before the cond mutex. */
    void (*parallel_routine)(int);
    std::vector<pthread_t> threads;
    
};

void *threadfunc(void *parm)
{
    synchronized_threads_t::thread_arg_t* args =
    (synchronized_threads_t::thread_arg_t*) parm;
    
    int tid = args->tid;
    
    pthread_mutex_lock(args->glob_mutex);
    if (args->init_routine) {
#ifdef _GNU_SOURCE
        VERBOSE(1,"Thread " << tid << " calling init routine on CPU "
                << (int) sched_getcpu() << endl);
#endif
        (*args->init_routine)(tid);
    }
    nthreads_initialized++;
    if (nthreads_initialized == nthreads) {
        wait_for_all_init = false;
        /* Release main thread */
        pthread_cond_signal(args->main_cond);
    }
    
    /* Wait for main to sync */
    while (!run_parallel) {
        pthread_cond_wait(args->cond, args->glob_mutex);
    }
#ifdef _GNU_SOURCE
    VERBOSE(1,"Thread " << tid << " calling parallel routine on CPU "
            << (int) sched_getcpu() << endl);
#endif
    
    pthread_mutex_unlock(args->glob_mutex);
    
    thread_timer[tid].start();
    
    if (args->parallel_routine) {
        (*args->parallel_routine)(tid);
    }
    
    thread_timer[tid].stop();
    
    return NULL;
}

void
synchronized_threads_t::execute_n_threads(unsigned n,
                                          void (*init_routine)(int),
                                          void (*parallel_routine)(int))
{
    nthreads = n;
    nthreads_initialized = 0;
    threads.resize(nthreads);
    thread_timer.resize(nthreads);
    
    for (int i=0; i<nthreads; i++) {
        synchronized_threads_t::thread_arg_t arg(i,init_routine,parallel_routine,
                                                 &glob_mutex, &cond, &main_cond);
        PZ_PTHREAD_CREATE(&threads[i],NULL,threadfunc,(void*) &i, __FUNCTION__);
    }
    
    /* Wait for all to be initialized */
    pthread_mutex_lock(&main_mutex);
    while (wait_for_all_init) {
        pthread_cond_wait(&main_cond, &main_mutex);
    }
    pthread_mutex_unlock(&main_mutex);
    
    /* Signall all to start together. */
    total_rst.start();
    run_parallel = true;
    pthread_cond_broadcast(&cond);
    
    /* Wait for all to finish. */
    for (unsigned i=0; i<nthreads; i++) {
        PZ_PTHREAD_JOIN(threads[i], NULL, __FUNCTION__);
    }
    total_rst.stop();
    
    if (verbose >= 2) {
        printf("%7s,%10s,%10s,%10s\n", "thread", "elapsed", "start", "stop");
        for (unsigned i=0; i<nthreads; i++) {
            printf("%7d,%10lld,%10lld,%10lld\n", i,
                   thread_timer[i].get_elapsed(),
                   thread_timer[i].get_start(),
                   thread_timer[i].get_stop());
        }
    }
}


int main(int argc, char *argv[])
{
#ifdef USING_TBB
    task_scheduler_init init;
#endif
    setup_masks();
    
    /* Parse the arguments */
    if (clarg::parse_arguments(argc, argv)) {
        cerr << "Error when parsing the arguments!" << endl;
        return 1;
    }
    
    verbose = verb_level.get_value();
    
    if (h.get_value() == true) {
        help(argv[0]);
        return 1;
    }
    
    if (nmats.get_value() < 1) {
        cerr << "Error, nmats must be >= 1" << endl;
        return 1;
    }
    
    if (verbose >= 1) {
        std::cout << "- Arguments -----------------------" << std::endl;
        clarg::values(std::cout, false);
        std::cout << "-----------------------------------" << std::endl;
    }
    
    synchronized_threads_t thread_exec;
    
    /* Read the matrix. */
    VERBOSE(1,"Reading input file: " << ifn.get_value() << std::endl);
    FileStreamWrapper input_file(bi.get_value());
    input_file.OpenRead(ifn.get_value());
    matrix.Read(input_file,0);
    VERBOSE(1,"Reading input file: " << ifn.get_value()
            << " [DONE]" << std::endl);
    
    if (maxcol.was_set())
        matrix.Resize(maxcol.get_value(),0);
    
    int nthreads = nmats.get_value();
    
    thread_exec.execute_n_threads(nthreads, init_decompose,
                                  compute_decompose);
    
    if (mstats.get_value() > 0) {
        unsigned n = matrix.Dim();
        uint64_t n_sky_items = 0;
        uint64_t max_height = 0;
        for (unsigned i=0; i<n; i++) {
            unsigned height = matrix.SkyHeight(i);
            if (mstats.get_value() > 1) {
                cout << "col " << i << " height = " << height << endl;
            }
            n_sky_items += height;
            if (height > max_height) max_height = height;
        }
        uint64_t n2 = n * n;
        double av_height = (double) n_sky_items / (double) n;
        cout << "N         = " << n << endl;
        cout << "N^2       = " << n2 << endl;
        cout << "Sky items = " << n_sky_items << endl;
        cout << "N^2 / Sky items = " << (double) n2 / (double) n_sky_items << endl;
        cout << "Avg. Height = " << av_height << endl;
        cout << "Max. Height = " << max_height << endl;
    }
    
    /** Dump decomposed matrix */
    if (dump_dm.was_set()) {
        VERBOSE(1, "Dumping decomposed matrix into: " <<
                dump_dm.get_value() << endl);
        FileStreamWrapper dump_file(bd.get_value());
        dump_file.OpenWrite(dump_dm.get_value());
        matrix.Write(dump_file, 0);
    }
    
    /* Gen/Check MD5 signature */
    if (gen_dm_sig.was_set() || chk_dm_sig.was_set()) {
        TPZMD5Stream sig;
        matrix.Write(sig, 1);
        int ret;
        if (chk_dm_sig.was_set()) {
            if ((ret=sig.CheckMD5(chk_dm_sig.get_value()))) {
                cerr << "ERROR(ret=" << ret << ") : MD5 Signature for "
                << "decomposed matrixdoes not match." << endl;
                return 1;
            }
            else {
                cout << "Checking decomposed matrix MD5 signature: [OK]" << endl;
            }
        }
        if (gen_dm_sig.was_set()) {
            if ((ret=sig.WriteMD5(gen_dm_sig.get_value()))) {
                cerr << "ERROR (ret=" << ret << ") when writing the "
                << "decomposed matrix MD5 signature to file: "
                << gen_dm_sig.get_value() << endl;
                return 1;
            }
        }
    }
    
    int ret=0; // Ok
    
    /** Check decomposed matrix */
    if (chk_dm_error.was_set()) {
        VERBOSE(1, "Checking decomposed matrix error: " <<
                chk_dm_error.get_value() << endl);
        FileStreamWrapper ref_file(br.get_value());
        ref_file.OpenRead(chk_dm_error.get_value());
        /* Reference matrix. */
        TPZSkylMatrix<REAL> ref_matrix;
        ref_matrix.Read(ref_file,0);
        int max_j = matrix.Cols();
        if (max_j != ref_matrix.Cols()) {
            cerr << "Decomposed matrix has " << max_j
            << " cols while reference matrix has "
            << ref_matrix.Cols() << endl;
            return 1;
        }
        REAL error_tolerance = error_tol.get_value();
        REAL max_error = 0.0;
        for (int j=0; j<max_j; j++) {
            int col_height = matrix.SkyHeight(j);
            if (col_height != ref_matrix.SkyHeight(j)) {
                cerr << "Column " << j << " of decomposed matrix has " << col_height
                << " non zero rows while reference matrix has "
                << ref_matrix.SkyHeight(j) << endl;
                return 1;
            }
            int min_i = (j+1) - col_height;
            for (int i=min_i; i<=j; i++) {
                
                REAL dm_ij = matrix.s(i,j);
                REAL rm_ij = ref_matrix.s(i,j);
                if (dm_ij != rm_ij) {
                    REAL diff = abs(dm_ij - rm_ij);
                    if (diff >= error_tolerance) {
                        VERBOSE(1, "diff(" << diff << ") tolerance (" << error_tolerance 
                                << "). dm[" << i << "][" << j << "] (" << dm_ij
                                << ") != rm[" << i << "][" << j << "] (" << rm_ij 
                                << ")." << endl);
                        ret = 1;
                        max_error = (max_error < diff)?diff:max_error;
                    }
                }
            }
        }
        if (ret != 0) {
            cerr << "Error ("<< max_error <<") > error tolerance ("
            << error_tolerance <<")" <<  endl;
        }
    }
    
    return ret;
}




