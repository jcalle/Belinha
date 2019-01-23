/**
 * @file
 * @brief Contains the implementation of the CG function which solves the symmetric positive definite linear system
 * using the Conjugate Gradient method.
 */

/**
 * @ingroup solver
 * @brief CG solves the symmetric positive definite linear system \f$ A x = b \f$ using the Conjugate Gradient method.
 * @return The return value indicates convergence within max_iter (input) iterations (0), or no convergence within max_iter iterations (1). \n
 * Upon successful return, output arguments have the following values: 
 * @param A  -- matrix of the system
 * @param b  -- vector of the system
 * @param M  -- preconditioner matrix
 * @param x  --  approximate solution to Ax = b
 * @param max_iter  --  the number of iterations performed before the tolerance was reached
 * @param tol  --  the residual after the final iteration
 * @param residual  -- residual vector (return)
 * @param FromCurrent  -- for type of operation (MultAdd)
*/
/**
 * Iterative template routine -- CG \n
 * CG follows the algorithm described on p. 15 in the SIAM Templates book.
 */
#ifdef TEST
#include <list> 
#endif

template < class Matrix, class Vector, class Preconditioner, class Real >
int
CG( Matrix &A, Vector &x, const Vector &b,
   Preconditioner &M, Vector *residual, int64_t &max_iter, Real &tol,const int FromCurrent)
{
	Real resid;
	Vector p, z, q;
	REAL alpha, beta, rho, rho_1 = 0;
	
    REAL normb = TPZExtractVal::val(Norm(b));
	Vector resbackup;
	Vector *res = residual;
	
#ifdef TEST
	std::list< TPZFMatrix<REAL> > plist,qlist;
	std::list< TPZFMatrix<REAL> >::iterator jt;
	std::list< TPZFMatrix<REAL> >::iterator kt;
	Vector Au;
#endif
	
	if(!res) res = &resbackup;
	Vector &r = *res;
	//  Vector r = b - A*x;
	if(FromCurrent) 
    {
        A.MultAdd(x,b,r,-1.,1.);
    }
	else {
		x.Zero();
		r = b;
	}
	
	if (normb == 0.0)
		normb = 1.0;
	
    if ((resid = (TPZExtractVal::val( Norm(r) ) ) / normb) <= tol) {
		tol = resid;
		max_iter = 0;
		return 0;
	}
	int64_t i;
	for (i = 1; i <= max_iter; i++) {
		M.Solve(r,z);
        rho = TPZExtractVal::val(Dot(r, z));
		
		if (i == 1)
			p = z;
		else {
			beta = rho / rho_1;
			p.TimesBetaPlusZ(beta,z);
		}
#ifdef TEST
		
		plist.push_back(p);
#endif	 
		
		A.Multiply(p,q);
		alpha = rho / (TPZExtractVal::val(Dot(p, q)));
		
#ifdef TEST
		qlist.push_back(q);
#endif
		
		x.ZAXPY(alpha,p);
		r.ZAXPY(-alpha,q);
		
#ifdef TEST
		A.Multiply(x,Au);
		REAL energy = Dot(x,Au)/2.-Dot(x,b);
#endif
		
		if ((resid = (TPZExtractVal::val(Norm(r))) / normb) <= tol) {
			tol = resid;
			max_iter = i;
#ifdef PZDEBUG
			std::cout << "cg iter = " << i <<  " res = " << resid << std::endl;
#endif
			return 0;
		}
#ifdef PZDEBUG
		std::cout << "cg iter = " << i <<  " res = " << resid /*<< " energy " << energy */ << std::endl;
#endif
#ifdef TEST
		std::cout << " energy " << energy << std::endl;
		TPZFMatrix<REAL> inner(plist.size(),plist.size(),0.);
		{
			int64_t j,k;
			for(j=0, jt = plist.begin(); jt != plist.end(); jt++,j++)
			{
				for(k=0, kt = qlist.begin(); kt != qlist.end(); kt++,k++)
				{
					inner(j,k) = Dot((*jt),(*kt));
				}
			}
		}
		inner.Print("Inner product of search directions");
#endif
		rho_1 = rho;
	}
	
	tol = resid;
	std::cout << "cg iter = " << i <<  " res = " << resid << std::endl;
	return 1;
}
