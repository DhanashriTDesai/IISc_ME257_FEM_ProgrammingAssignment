#include <iostream>
#include <cmath>
#include <Eigen/Sparse> 
#include <fstream>
#include <Vector>
#include "ReadMesh.h"
#include "WriteMesh.h"


using namespace std;
using namespace Eigen;


// Declare and define - Shape functions and derivatives over the PARAMETRIC element: 
// xi => Input i.e Coordinates in the parametric triangle 
// a => Input i.e. Shape function number, possible values 0, 1, 2
// val => Output i.e. Value of the shape function 
// dval => Output i.e. Derivatives of the shape function 

void GetParametricShapeFunction(double* xi, int a, double& valpara, double* dvalpara)
{
	if (a == 0)
	{
		valpara = 1 - xi[0] - xi[1];
		dvalpara[0] = -1;
		dvalpara[1] = -1;
	}

	else if (a == 1)
	{
		valpara = xi[0];
		dvalpara[0] = 1;
		dvalpara[1] = 0;
	}

	else if (a == 2)
	{
		valpara = xi[1];
		dvalpara[0] = 0;
		dvalpara[1] = 1;
	}

	// Print Value of Shape function
	//cout << "\nValue of the Parametric Shape Function";
	//cout << "\n" << valpara << "\n";

	// Print Value of Gradient of Shape function
	//cout << "\nValue of Gradient of Parametric Shape Function";
	//for (int i = 0; i < 2 ; i++)
		//cout << "\n" << dvalpara[i];
	//cout << "\n";
}


// Declare and define - Shape functions and derivatives over the PHYSICAL element: 
// xi => Input i.e. Coordinates of a quadrature point in the parametric triang e
// a => Input i.e. Local shape function number in an element Possible values : 0, 1, 2
// nodal_coords => Input i.e. Cartesian coordinates of the 3 nodes of a triangle element
// val => Output i.e. Value of the shape function Na
// dval => Output i.e. Derivative of the shape function Na

void GetPhysicalShapeFunction(double* xi, int a, double* nodal_coords, double& val, double* dval, double e, double& det)
{
	double dvalpara[2];
	double valpara = 0;

	double del_phi[2][2] = { 0,0,0,0 };
	del_phi[0][0] = nodal_coords[2] - nodal_coords[0];
	del_phi[1][0] = nodal_coords[3] - nodal_coords[1];
	del_phi[0][1] = nodal_coords[4] - nodal_coords[0];
	del_phi[1][1] = nodal_coords[5] - nodal_coords[1];

	double del_phiT[2][2] = { 0,0,0,0 };
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			del_phiT[i][j] = del_phi[j][i];

	det = del_phiT[0][0] * del_phiT[1][1] - del_phiT[1][0] * del_phiT[0][1];

	double Adjoint[2][2] = { 0,0,0,0 };
	Adjoint[0][0] = del_phiT[1][1];
	Adjoint[1][1] = del_phiT[0][0];
	Adjoint[0][1] = -del_phiT[0][1];
	Adjoint[1][0] = -del_phiT[1][0];

	double del_phiTinv[2][2] = { 0,0,0,0 };
	del_phiTinv[0][0] = (1 / det) * Adjoint[0][0];
	del_phiTinv[1][1] = (1 / det) * Adjoint[1][1];
	del_phiTinv[0][1] = (1 / det) * Adjoint[0][1];
	del_phiTinv[1][0] = (1 / det) * Adjoint[1][0];

	//cout << "\n\nJacobian for Element " << e << " is : \n";
	//for (int i = 0; i < 2; i++)
	//{
		//for (int j = 0; j < 2; j++)
			//cout << del_phi[i][j] << "\t";
		//cout << "\n";
	//}

	//cout << "\nTranspose of Jacobian for Element " << e << " is : \n";
	//for (int i = 0; i < 2; i++)
	//{
		//for (int j = 0; j < 2; j++)
			//cout << del_phiT[i][j] << "\t";
		//cout << "\n";
	//}

	//cout << "\n\nInverse Transpose of Jacobian for Element " << e << " is : \n";
	//for (int i = 0; i < 2; i++)
	//{
		//for (int j = 0; j < 2; j++)
			//cout << del_phiTinv[i][j] << "\t";
		//cout << "\n";
	//}

	GetParametricShapeFunction(xi, a, valpara, dvalpara);

	val = valpara;
	//cout << "\nValue of Shape Function is " << val;

	//cout << "\n\nGradient of Parametric Shape function for node : " << a;
	//cout << "\ndvalpara[0] : " << dvalpara[0] << "\ndvalpara[1] : " << dvalpara[1];

	dval[0] = del_phiTinv[0][0] * dvalpara[0] + del_phiTinv[0][1] * dvalpara[1];
	dval[1] = del_phiTinv[1][0] * dvalpara[0] + del_phiTinv[1][1] * dvalpara[1];
	//cout << "\n\nGradient of Physical Shape function for node : " << a;
	//cout << "\ndval[0] : " << dval[0] << "\ndval[1] : " << dval[1];
}


// Function to compute 3x3 Element Stiffness Matrix
// nodal coords => Input i.e. Nodal coordinates of the element
// kmat => Output Element i.e. stiffness matrix Size 3x3

void GetElementStiffnessMatrix(double* nodal_coords, MatrixXd& ke, double e, double det)
{

	double w = 0.1667;
	double quadpts[] = { 0.5 , 0 , 0.5 , 0.5 , 0 , 0.5 };

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int q = 0; q < 3; q++)
			{
				double val;
				double dvalone[2] = { 0,0 };
				double dvaltwo[2] = { 0,0 };
				double xi[] = { quadpts[2 * q] ,quadpts[2 * q + 1] };

				GetPhysicalShapeFunction(xi, i, nodal_coords, val, dvalone, e, det);
				GetPhysicalShapeFunction(xi, j, nodal_coords, val, dvaltwo, e, det);

				ke(i, j) += w * det * (dvalone[0] * dvaltwo[0] + dvalone[1] * dvaltwo[1]);
			}
		}
	}

	//cout << "\n\n Element Stffness Matrix for Element " << e;
	//cout << "\n" << ke << "\n";
}


// Function to compute 3x1 Element force vector
// nodal coords => Input i.e. Nodal coordinates of the element
// fvec => Output i.e Element force vector

void GetElementForceVector(double* nodal_coords, double* fvec, double e, double det)
{
	double w = 0.1667;
	double quadpts[] = { 0.5 , 0 , 0.5 , 0.5 , 0 , 0.5 };

	for (int i = 0; i < 3; i++)
	{
		for (int q = 0; q < 3; q++)
		{
			double val;
			double dval[2] = { 0,0 };
			double xi[] = { quadpts[2 * q] ,quadpts[2 * q + 1] };
			double x, y, f;
			const double pi = 3.142;

			GetPhysicalShapeFunction(xi, i, nodal_coords, val, dval, e, det);
			x = nodal_coords[0] * (1 - xi[0] - xi[1]) + nodal_coords[2] * (xi[0]) + nodal_coords[4] * (xi[1]);
			y = nodal_coords[1] * (1 - xi[0] - xi[1]) + nodal_coords[3] * (xi[0]) + nodal_coords[5] * (xi[1]);
			f = 0;
			fvec[i] += w * det * val * f;
		}
	}

	//cout << "\n Element Force Vector for Element " << e << "\n";
	//cout << "\n" << fvec[0] << "\n" << fvec[1] << "\n" << fvec[2] << "\n";
}


// Funtion to  Write Important things to file - 'Solution-HW5-Trial.txt'
void WriteToFile(int nNodes , int nElements , vector<double> coordinates , vector<int> connectivity , SparseMatrix<double, RowMajor> K , VectorXd F, VectorXd alpha)
{
	ofstream outFile;
	outFile.open("Solution-HW5-Trial.txt");
	outFile << "Read " << nNodes << " nodes and " << nElements << " elements ";

	// Write Coordinates Array to file
	outFile << "\n\nCoordinates Array: \n";
	for (int i = 0; i < coordinates.size(); i++)
		outFile << coordinates[i] << "  ";

	// Write Connectivity Array to file
	outFile << "\n\nConnectivity Array: \n";
	for (int i = 0; i < connectivity.size(); i++)
		outFile << connectivity[i] << "  ";

	// Write Global Stiffness matrix to file
	outFile << "\n\nGlobal Stiffness matrix is : ";
	outFile << "\n" << MatrixXd(K) << "\n";

	// Write Global Force Vector to file
	outFile << "\n\nGlobal Force Vector is : ";
	outFile << "\n" << F << "\n";

	// Write Solution Vector to file
	outFile << "\n\nSolution Vector is : ";
	outFile << "\n" << alpha << "\n";
	outFile.close();
}


int main()
{
	// Read-in the mesh
	vector<double> coordinates;
	vector<int> connectivity;
	ReadMesh("coordinates-5.dat", coordinates, "connectivity-5.dat", connectivity);		// Reads in the Mesh
	PlotTecFile("VisualizeMesh.tec", coordinates, connectivity);		// Creates PlotTecFile - 'VisualizeMesh.tec' file for Mesh Visualization

	// Compute nNodes and nElements from sizes of coordinates and connectivity array respectively
	int nNodes = coordinates.size() / 2;
	int nElements = connectivity.size() / 3;

	int a;
	double det = 0;
	double nodal_coords[6] = { 0,0,0,0,0,0 };

	// Define Vector of type VectorXd for the Global Force Vector 'F'
	VectorXd F(nNodes);			
	F = VectorXd::Zero(nNodes);		// Initialize to zero

	// Define Matrix of type SparseMatrix<double, RowMajor> for the Global Stiffness matrix 'K'
	SparseMatrix<double, RowMajor> K(nNodes, nNodes);
	K.setZero();		// Initialize to zero

	// Define a Vector 'tripletVector' to store NON-ZERO entries of Sparse Matrix 'gMassMatrix'
	typedef Triplet<double> T;
	vector<T> tripletVector;

	for (int e = 0; e < nElements; e++)
	{
		for (a = 0; a < 3; a++)
		{
			int l2gmap;
			l2gmap = connectivity[3 * e + a];

			nodal_coords[2 * a + 0] = coordinates[2 * l2gmap + 0];
			nodal_coords[2 * a + 1] = coordinates[2 * l2gmap + 1];
		}

		// Define Matrix of type MatrixXd for the Element Stiffness matrix 'ke'
		MatrixXd ke(3, 3);	
		ke.setZero();		// Initialize to zero

		GetElementStiffnessMatrix(nodal_coords, ke, e, det);

		int p = connectivity[3 * e];
		int q = connectivity[3 * e + 1];
		int r = connectivity[3 * e + 2];

		tripletVector.push_back(T(p, p, ke(0, 0)));
		tripletVector.push_back(T(p, q, ke(0, 1)));
		tripletVector.push_back(T(p, r, ke(0, 2)));
		tripletVector.push_back(T(q, p, ke(1, 0)));
		tripletVector.push_back(T(q, q, ke(1, 1)));
		tripletVector.push_back(T(q, r, ke(1, 2)));
		tripletVector.push_back(T(r, p, ke(2, 0)));
		tripletVector.push_back(T(r, q, ke(2, 1)));
		tripletVector.push_back(T(r, r, ke(2, 2)));

		// Define Element Force Vector 'fvec' and Initialize to Zero
		double fvec[3] = { 0,0,0 };

		GetElementForceVector(nodal_coords, fvec, e, det);

		// Assemble Global Force Vector using Element Force Vector 'fvec'
		F(p) += fvec[0];
		F(q) += fvec[1];
		F(r) += fvec[2];
	}

	// Assemble Global Stiffness Matrix 'K' using 'setFromTriplets()' method.
	K.setFromTriplets(tripletVector.begin(), tripletVector.end());

	// Enforce Non Homo BCs (u=y along boundary) on Global force Vector F and Global Stiffness Matrix K
	for (int k = 0; k < coordinates.size() / 2; k++)
	{
		if (coordinates[2 * k] == -1 || coordinates[2 * k] == 1 || coordinates[2 * k + 1] == -1 || coordinates[2 * k + 1] == 1)
		{
			F(k) = coordinates[2 * k + 1];		// ensures u=x along boundary
			for (int j = 0; j < nNodes; j++)
				K.coeffRef(k, j) = 0;
			K.coeffRef(k, k) = 1;
		}
	}

	// Solve Linear System K(alpha)=F
	VectorXd alpha(nNodes);		// Define Solution Vector 'alpha'
	alpha = VectorXd::Zero(nNodes);		// Initialize to zero

	// Use SparseLU Routine to solve Linear System
	SparseLU<SparseMatrix<double>> solver;
	solver.analyzePattern(K);
	solver.factorize(K);
	alpha = solver.solve(F);

	// Create 1D Array 'solution' for Generating PlotTecFile file of Solution to visualise FE Solution
	double solution[10000];
	for (int i = 0; i < nNodes; i++)
		solution[i] = alpha(i);

	// Create PlotTecFile - 'VisualizeSolution.tec' file of 'solution' 
	PlotTecFile("VisualizeSolution.tec", coordinates, connectivity, solution);		

	WriteToFile(nNodes, nElements, coordinates, connectivity, K, F, alpha);

	return 0;
}

