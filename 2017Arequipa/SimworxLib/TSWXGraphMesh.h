//---------------------------------------------------------------------------

#ifndef TSWXGraphMeshH
#define TSWXGraphMeshH

//#include "TSWXVector.h"
#include "pzvec.h"
#include "pzstack.h"
#include "string"
#include <map>

enum EGraphElType{ ENoneGraphElType = 0,
                   EvtkQuadraticEdge = 10,
                   EvtkBiQuadraticTriangle = 20,
                   EvtkBiQuadraticQuad = 30,
                   EvtkQuadraticTetra = 40,
                   EvtkQuadraticPyramid = 50,
                   EvtkQuadraticWedge = 60,
                   EvtkQuadraticHexahedron = 70,
                   EvtkLinearHexahedron = 80
									 };

enum ESWXVtkSolType{ ECellSolution = 1, ENodeSolution = 2 };

typedef TPZVec< double > TSWXGraphNode;

typedef TPZVec< int > TSWXGraphIncid;

struct TSWXGraphEl{

public:
  TSWXGraphIncid fIncid;
  EGraphElType fElType;
  int fMatId;

  void Read(std::istream & file);
  void Write(std::ostream &file) const;

};

typedef TPZVec<double> TSWXGraphSingleSol;

struct TSWXGraphSol{

public:
  ESWXVtkSolType fSolType;
  std::string fTitle;
  TPZStack< TSWXGraphSingleSol > fData;
	int Dimension() const;
	unsigned int Size() const;
  void Read(std::istream & file);
  void Write(std::ostream &file) const;
};

/** Guarda solu��o para diversos passos de tempo
 */
class TSWXGraphMeshSol{

private:

	///Mapa de tempo para solu��o
	std::map< double, TPZVec< TSWXGraphSol > > fAllSol;

public:

	///Construtor
	TSWXGraphMeshSol();

	///Destrutor
	virtual ~TSWXGraphMeshSol();

	///Adiciona solu��o para um passo de tempo
	void AddSol(double time, const TPZVec< TSWXGraphSol >  &sol);

	///Retorna o n�mero de passos de tempo dos resultados
	unsigned int NTimeSteps() const;

	///Retorna a solu��o de um passo de tempo
	void GetSolData(int iStep, double &time, TPZVec< TSWXGraphSol >  &sol) const;

	///Retorna a solu��o de um passo de tempo
	const TPZVec< TSWXGraphSol >  & GetSolData(int iStep, double &time) const;

	///Apaga todos os resultados salvos
	void Clear();

	TPZVec<std::string> GetSolutionsNames() const;
	void GetSolution(std::string SolName, std::map<double , TSWXGraphSol> & solMap) const;

	void Read(std::istream & file);

	void Write(std::ostream &file) const;

};

class TSWXGraphMesh{

public:
	TPZStack< TSWXGraphNode > fNodes;
	TPZStack< TSWXGraphEl > fElem;
	TSWXGraphMeshSol fSol;
  std::map< int, std::string > fMaterialLabels;

  void Read(std::istream & file);
  void Write(std::ostream &file) const;

//  _di_IXMLNode ReadMe(_di_IXMLNode & fatherNode, const UnicodeString & nodeName);
//  _di_IXMLNode WriteMe(_di_IXMLNode & fatherNode, const UnicodeString & nodeName) const;

  ///Reset all data
  void Reset(){
    fNodes.resize(0);
    fElem.resize(0);
    fSol.Clear();
    fMaterialLabels.clear();
	}

	bool HasData() const{ return (fElem.size()); }

  /** Escreve arquivo paraview da solu��o do passo de tempo
   * mathStyleIStep. Se mathStyleIStep >= 0, ent�o � escrita
   * a solu��o do passo de tempo iStep.
   * Se mathStyleIStep < 0, ent�o � escrita a solu��o do passo
   * de tempo fSol.NTimeSteps() - mathStyleIStep, i.e.,
   * se mathStyleIStep == -1 a �ltima solu��o � escrita.
   * Foi feito assim para manter os c�digos j� existentes.
   */
  void ToParaview(std::ostream &file, int mathStyleIStep = -1) const;

  bool HasNodalSolution() const;
  bool HasCellSolution() const;
  int VTKCellType( EGraphElType enumType ) const;
};

#endif
