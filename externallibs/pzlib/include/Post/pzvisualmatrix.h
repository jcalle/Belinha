/**
 * @file
 * @brief Contains the declaration of the VisualMatrix functions to VTK and DX packages.
 */

#ifndef VISUAL_MATRIX
#define VISUAL_MATRIX

#include "pzfmatrix.h"
#include "pzvec.h"
#include "fstream"

/** 
 * @ingroup post
 * @brief This function calls the function that create a Data Explorer file or \n
 * VTK file that allow to visualization of the value of a matrix passed as parameter.
 */
template<class TVar>
void VisualMatrix(TPZFMatrix<TVar> &matrix, const std::string &outfilename);

/** 
 * @ingroup post
 * @brief This function creates a Data Explorer file that allow to visualization of the value of a matrix passed as parameter.
 */
template<class TVar>
void VisualMatrixDX(TPZFMatrix<TVar> &matrix, const std::string &outfilename);

/** 
 * @ingroup post
 * @brief This function creates a VTK file that allow to visualization of the value of a matrix passed as parameter.
 */
template<class TVar>
void VisualMatrixVTK(TPZFMatrix<TVar> &matrix, const std::string &outfilename);

#endif
