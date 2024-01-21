#ifndef TRIANGLEREFINER_HPP
#define TRIANGLEREFINER_HPP

#include "GenericDomain.hpp"
#include "GenericMesh.hpp"
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace GeDiM
{
	class TriangleRefiner
	{
	protected:
		GenericMesh* meshPointer;
		vector<unsigned int> idCellsToRefine;
	public:
		TriangleRefiner() { meshPointer = NULL; }
		~TriangleRefiner() { }

		void SetMesh( GenericMesh& mesh ) { meshPointer = &mesh; }
		void SetNumberCellsToRefine( const unsigned int& value ) { idCellsToRefine.reserve(value); }
		void AddCellId( const unsigned int& value ) { idCellsToRefine.push_back(value);}
		void Rotate( const unsigned int& value );
		void UpdateNeighbourhood();
		void PrintNeigh();



	  Output::ExitCodes RefineMesh();
	  Output::ExitCodes RefineTriangle( const unsigned int& value ); //Riceve l'indice del triangolo da raffinare
	  Output::ExitCodes RefineUniformly( const unsigned int& value ); //Taglia il triangolo in 4 triangoli simili
	};
}

#endif
