#ifndef MESHIMPORT_TRIANGLE_H
#define MESHIMPORT_TRIANGLE_H

#include <vector>

#include "Eigen/Eigen"

#include "GenericMesh.hpp"
#include "GenericDomain.hpp"

#include "mytriangle.hpp"

using namespace std;
using namespace Eigen;
using namespace MainApplication;

namespace GeDiM
{
	class MeshImport_Triangle;

	class MeshImport_Triangle : public GenericMeshImportInterface
	{
		protected:
			struct triangulateio* inputMeshPointer; //usa le struct perchè scritto originariamente in C
			struct triangulateio* outputMeshPointer;

			string triangleOptions;

		public:
			MeshImport_Triangle();
			virtual ~MeshImport_Triangle();

			void SetInputMesh(struct triangulateio* _inputMeshPointer) { inputMeshPointer = _inputMeshPointer; }
			void SetTriangleOptions(const string& _triangleOptions) { triangleOptions = _triangleOptions; }

			static int CheckTrianglePosition(const vector<int>& edgePointIds,const vector<int>& cellPointIds);

			Output::ExitCodes CreateTriangleInput(const GenericDomain& domain);
			Output::ExitCodes CreateTriangleOutput(const GenericDomain& domain);
			Output::ExitCodes CreateMesh(const GenericDomain& domain, GenericMesh& mesh) const;

			Output::ExitCodes ExportTriangleMesh(const string& nameFolder, const string& nameFile) const;
	};
}

#endif // MESHIMPORT_TRIANGLE_H
