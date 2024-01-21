#ifndef GENERICMESH_H
#define GENERICMESH_H

#include <vector>
#include <map>

#include "Eigen/Eigen"
#include "Output.hpp"
#include "GenericDomain.hpp"

using namespace std;
using namespace Eigen;
using namespace MainApplication;

namespace GeDiM
{
	class GenericDomain;
	class GenericTreeNode;
	class GenericCell;
	class GenericFace;
	class GenericEdge;
	class GenericPoint;
	class GenericMesh;
	class GenericMeshImportInterface;

	class GenericTreeNode
	{
		protected:
			const GenericTreeNode* father; ///< Father node in the tree
			vector<const GenericTreeNode*> childs; ///< Child nodes in the tree.

			unsigned int id; ///< Internal id of the node
			unsigned int globalId; ///< Global id of the node, generally equal to id
			bool isActive; ///< Tells if the node is active for system matrix computation

			unsigned short marker; ///< Marker for Dirichlet and Neumann conditions

			map<string, void*> properties;

		public:
			GenericTreeNode(const unsigned int& _id);
			GenericTreeNode(const GenericTreeNode& treeNode);
			virtual ~GenericTreeNode();

			const unsigned int& Id() const { return id; }
			const unsigned int& GlobalId() const { return globalId; }
			const bool& IsActive() const { return isActive; }
			const unsigned short& Marker() const { return marker; }
			const bool IsDirichlet() const { return (marker % 2 == 1); }
			const bool IsNeumann() const { return (marker > 0 && marker % 2 == 0); }
			const bool HasFather() const { return (father != NULL);}
			const bool HasChilds() const { return (childs.size() > 0);}
			const unsigned int NumberOfChilds() const { return childs.size();}

			const GenericTreeNode* Father() const{ return father;}
			const GenericTreeNode* Child(const unsigned int& position) const{ return childs[position];}

			void SetGlobalId(const unsigned int& _globalId) { globalId = _globalId; }
			void SetId(const unsigned int& _id) { id = _id; }
			void SetState(bool _isActive = true) { isActive = _isActive; }
			void SetMarker(const unsigned short& _marker) { marker = _marker; }

			void SetFather(const GenericTreeNode* _father) {father = _father;}

			void AllocateChilds(const unsigned int& numChilds){ childs.resize( numChilds);}
			void InsertChild(const GenericTreeNode* child, unsigned int& position) {childs[position] = child;}

			void InitializeChilds(const unsigned int& numChilds){ childs.reserve( numChilds);}
			void AddChild(const GenericTreeNode* child) {childs.push_back(child);}
			Output::ExitCodes InheritPropertiesByFather();

			const Output::ExitCodes InitializeProperty(const string& key);
			void AddProperty(const string& key, void* value) { properties[key] = value;}
			const map<string, void* >& GetAllProperties() const { return properties; }
			map<string, void* >& GetAllProperties() { return properties; }
			const void* GetProperty(const string& key) const { return properties.at(key); }
			void* GetProperty(const string& key) { return properties.at(key); }
			const bool HasProperty(const string& key) const { return (properties.find(key) != properties.end()); }
	};

	class GenericCell : public GenericTreeNode
	{
		protected:
			vector<const GenericCell*> cells; ///< Cells of the cell
			vector<const GenericFace*> faces; ///< Faces of the cell
			vector<const GenericEdge*> edges; ///< Edges of the cell
			vector<const GenericPoint*> points; ///< Points of the cell

		public:
			GenericCell(const unsigned int& _id);
			GenericCell(const GenericCell& cell);
			virtual ~GenericCell();

			const size_t NumberOfCells() const { return cells.size(); }
			const size_t NumberOfFaces() const { return faces.size(); }
			const size_t NumberOfEdges() const { return edges.size(); }
			const size_t NumberOfPoints() const { return points.size(); }
			const GenericCell* Cell(const unsigned int& position) const { return cells[position]; }
			const GenericFace* Face(const unsigned int& position) const { return faces[position]; }
			const GenericEdge* Edge(const unsigned int& position) const { return edges[position]; }
			const GenericPoint* Point(const unsigned int& position) const { return points[position]; }

			void InitializeCells(const size_t numberOfCells) { cells.reserve(numberOfCells); }
			void InitializeFaces(const size_t numberOfFaces) { faces.reserve(numberOfFaces); }
			void InitializeEdges(const size_t numberOfEdges) { edges.reserve(numberOfEdges); }
			void InitializePoints(const size_t numberOfPoints) { points.reserve(numberOfPoints); }

			void AllocateCells(const size_t numberOfCells) { cells.resize(numberOfCells, NULL); }
			void AllocateFaces(const size_t numberOfFaces) { faces.resize(numberOfFaces, NULL); }
			void AllocateEdges(const size_t numberOfEdges) { edges.resize(numberOfEdges, NULL); }
			void AllocatePoints(const size_t numberOfPoints) { points.resize(numberOfPoints, NULL); }

			void ErasePoint(const unsigned int& position) { points.erase(points.begin() + position); }
			void EraseEdge(const unsigned int& position) { edges.erase(edges.begin() + position); }
			void EraseFace(const unsigned int& position) { faces.erase(faces.begin() +position); }
			void EraseCell(const unsigned int& position) { cells.erase(cells.begin() + position); }

			void ShrinkCells() { cells.shrink_to_fit(); }
			void ShrinkFaces() { faces.shrink_to_fit(); }
			void ShrinkEdges() { edges.shrink_to_fit(); }
			void ShrinkPoints() { points.shrink_to_fit(); }

			Output::ExitCodes AddCell(const GenericCell* cell);
			Output::ExitCodes AddFace(const GenericFace* face);
			Output::ExitCodes AddEdge(const GenericEdge* edge);
			Output::ExitCodes AddPoint(const GenericPoint* point);

			Output::ExitCodes InsertCell(const GenericCell* cell, const unsigned int& position);
			Output::ExitCodes InsertFace(const GenericFace* face, const unsigned int& position);
			Output::ExitCodes InsertEdge(const GenericEdge* edge, const unsigned int& position);
			Output::ExitCodes InsertPoint(const GenericPoint* point, const unsigned int& position);

			bool PointInCellAndIdEdgeBoundary (const Vector3d& point, int& numEdges, const double& toll  = 1.0E-7) const;
			bool PointInCell(const Vector3d& point, const double& toll = 1.0E-7) const ;
	};

	class GenericFace : public GenericTreeNode
	{
		protected:
			vector<const GenericCell*> cells; ///< Cells of the face
			vector<const GenericFace*> faces; ///< Faces of the face
			vector<const GenericEdge*> edges; ///< Edges of the face
			vector<const GenericPoint*> points; ///< Points of the face

		public:
			GenericFace(const unsigned int& _id);
			GenericFace(const GenericFace& face);
			virtual ~GenericFace();

			const size_t NumberOfCells() const { return cells.size(); }
			const size_t NumberOfFaces() const { return faces.size(); }
			const size_t NumberOfEdges() const { return edges.size(); }
			const size_t NumberOfPoints() const { return points.size(); }
			const GenericCell* Cell(const unsigned int& position) const { return cells[position]; }
			const GenericFace* Face(const unsigned int& position) const { return faces[position]; }
			const GenericEdge* Edge(const unsigned int& position) const { return edges[position]; }
			const GenericPoint* Point(const unsigned int& position) const { return points[position]; }

			void InitializeCells(const size_t numberOfCells) { cells.reserve(numberOfCells); }
			void InitializeFaces(const size_t numberOfFaces) { faces.reserve(numberOfFaces); }
			void InitializeEdges(const size_t numberOfEdges) { edges.reserve(numberOfEdges); }
			void InitializePoints(const size_t numberOfPoints) { points.reserve(numberOfPoints); }

			void AllocateCells(const size_t numberOfCells) { cells.resize(numberOfCells, NULL); }
			void AllocateFaces(const size_t numberOfFaces) { faces.resize(numberOfFaces, NULL); }
			void AllocateEdges(const size_t numberOfEdges) { edges.resize(numberOfEdges, NULL); }
			void AllocatePoints(const size_t numberOfPoints) { points.resize(numberOfPoints, NULL); }

			void ErasePoint(const unsigned int& position) {points.erase(points.begin() + position);}
			void EraseEdge(const unsigned int& position) {edges.erase(edges.begin() + position);}
			void EraseFace(const unsigned int& position) {faces.erase(faces.begin() +position);}
			void EraseCell(const unsigned int& position) {cells.erase(cells.begin() + position);}

			void ShrinkCells() { cells.shrink_to_fit(); }
			void ShrinkFaces() { faces.shrink_to_fit(); }
			void ShrinkEdges() { edges.shrink_to_fit(); }
			void ShrinkPoints() { points.shrink_to_fit(); }

			Output::ExitCodes AddCell(const GenericCell* cell);
			Output::ExitCodes AddFace(const GenericFace* face);
			Output::ExitCodes AddEdge(const GenericEdge* edge);
			Output::ExitCodes AddPoint(const GenericPoint* point);

			Output::ExitCodes InsertCell(const GenericCell* cell, const unsigned int& position);
			Output::ExitCodes InsertFace(const GenericFace* face, const unsigned int& position);
			Output::ExitCodes InsertEdge(const GenericEdge* edge, const unsigned int& position);
			Output::ExitCodes InsertPoint(const GenericPoint* point, const unsigned int& position);

			const Matrix3d RotationMatrix(double rotationTolerance = 1.0E-12);
	};

	class GenericEdge : public GenericTreeNode
	{
		public:
			enum PositionPoint
			{
				AtTheLeft = 0,
				AtTheRight = 1,
				Beyond = 2,
				Behind = 3,
				Between = 4,
				AtTheOrigin = 5,
				AtTheEnd = 6
			};

		protected:
			vector<const GenericCell*> cells; ///< Cells of the edge
			vector<const GenericFace*> faces; ///< Faces of the edge
			vector<const GenericEdge*> edges; ///< Edges of the edge
			vector<const GenericPoint*> points; ///< Points of the edge

		public:
			GenericEdge(const unsigned int& _id);
			GenericEdge(const unsigned int& _id, const GenericPoint* origin, const GenericPoint* end);
			GenericEdge(const GenericEdge& edge);
			virtual ~GenericEdge();

			const size_t NumberOfCells() const { return cells.size(); }
			const size_t NumberOfFaces() const { return faces.size(); }
			const size_t NumberOfEdges() const { return edges.size(); }
			const size_t NumberOfPoints() const { return points.size(); }

			//position 0 = right position 1 = left
			const GenericCell* Cell(const unsigned int& position) const { return cells[position]; }
			const GenericCell* RightCell() const {return cells[0];}
			const GenericCell* LeftCell() const {return cells[1];}
			const GenericFace* Face(const unsigned int& position) const { return faces[position]; }
			const GenericEdge* Edge(const unsigned int& position) const { return edges[position]; }
			const GenericPoint* Point(const unsigned int& position) const { return points[position]; }

			const bool HasRightCell() const { return (cells[0] != NULL); }
			const bool HasLeftCell() const { return (cells[1] != NULL); }

			void InitializeCells(const size_t numberOfCells) { cells.reserve(numberOfCells); }
			void InitializeFaces(const size_t numberOfFaces) { faces.reserve(numberOfFaces); }
			void InitializeEdges(const size_t numberOfEdges) { edges.reserve(numberOfEdges); }

			void AllocateCells(const size_t numberOfCells) { cells.resize(numberOfCells, NULL); }
			void AllocateFaces(const size_t numberOfFaces) { faces.resize(numberOfFaces, NULL); }
			void AllocateEdges(const size_t numberOfEdges) { edges.resize(numberOfEdges, NULL); }
			void AllocatePoints(const size_t numberOfPoints) { points.resize(numberOfPoints, NULL); }

			void ErasePoint(const unsigned int& position) {points.erase(points.begin() + position);}
			void EraseEdge(const unsigned int& position) {edges.erase(edges.begin() + position);}
			void EraseFace(const unsigned int& position) {faces.erase(faces.begin() +position);}
			void EraseCell(const unsigned int& position) {cells.erase(cells.begin() + position);}

			void ShrinkCells() { cells.shrink_to_fit(); }
			void ShrinkFaces() { faces.shrink_to_fit(); }
			void ShrinkEdges() { edges.shrink_to_fit(); }
			void ShrinkPoints() { points.shrink_to_fit(); }

			Output::ExitCodes AddCell(const GenericCell* cell);
			Output::ExitCodes AddFace(const GenericFace* face);
			Output::ExitCodes AddEdge(const GenericEdge* edge);
			Output::ExitCodes AddPoint(const GenericPoint* point);

			Output::ExitCodes InsertCell(const GenericCell* cell, const unsigned int& position);
			Output::ExitCodes InsertFace(const GenericFace* face, const unsigned int& position);
			Output::ExitCodes InsertEdge(const GenericEdge* edge, const unsigned int& position);
			Output::ExitCodes InsertPoint(const GenericPoint* point, const unsigned int& position);

			Output::ExitCodes ChangeOrientation();
			GenericEdge::PositionPoint PointOnEdge (const Vector3d& point, const double& toll = 1.0E-7) const;
	};

	class GenericPoint : public GenericTreeNode
	{
		protected:
			Vector3d coordinates; ///< Geometry vertex of the degree of freedom. Its size depends on the dimension

			vector<const GenericCell*> cells; ///< Cells of the point
			vector<const GenericFace*> faces; ///< Faces of the point
			vector<const GenericEdge*> edges; ///< Edges of the point

		public:
			GenericPoint(const unsigned int& _id);
			GenericPoint(const GenericPoint& point);
			virtual ~GenericPoint();

			Vector3d& Coordinates() { return coordinates; }
			const Vector3d& Coordinates() const { return coordinates; }
			const double& X() const { return coordinates[0]; }
			const double& Y() const { return coordinates[1]; }
			const double& Z() const { return coordinates[2]; }
			const size_t NumberOfCells() const { return cells.size(); }
			const size_t NumberOfFaces() const { return faces.size(); }
			const size_t NumberOfEdges() const { return edges.size(); }
			const GenericCell* Cell(const unsigned int& position) const { return cells[position]; }
			const GenericFace* Face(const unsigned int& position) const { return faces[position]; }
			const GenericEdge* Edge(const unsigned int& position) const { return edges[position]; }

			Output::ExitCodes SetCoordinates(const Vector3d& _coordinates);
			Output::ExitCodes SetCoordinates(const double& x, const double& y = 0.0, const double& z = 0.0);

			void InitializeCells(const size_t numberOfCells) { cells.reserve(numberOfCells); }
			void InitializeFaces(const size_t numberOfFaces) { faces.reserve(numberOfFaces); }
			void InitializeEdges(const size_t numberOfEdges) { edges.reserve(numberOfEdges); }

			void AllocateCells(const size_t numberOfCells) { cells.resize(numberOfCells, NULL); }
			void AllocateFaces(const size_t numberOfFaces) { faces.resize(numberOfFaces, NULL); }
			void AllocateEdges(const size_t numberOfEdges) { edges.resize(numberOfEdges, NULL); }

			void EraseEdge(const unsigned int& position) {edges.erase(edges.begin() + position);}
			void EraseFace(const unsigned int& position) {faces.erase(faces.begin() +position);}
			void EraseCell(const unsigned int& position) {cells.erase(cells.begin() + position);}

			void ShrinkCells() { cells.shrink_to_fit(); }
			void ShrinkFaces() { faces.shrink_to_fit(); }
			void ShrinkEdges() { edges.shrink_to_fit(); }

			Output::ExitCodes AddCell(const GenericCell* cell);
			Output::ExitCodes AddFace(const GenericFace* face);
			Output::ExitCodes AddEdge(const GenericEdge* edge);

			Output::ExitCodes InsertCell(const GenericCell* cell, const unsigned int& position);
			Output::ExitCodes InsertFace(const GenericFace* face, const unsigned int& position);
			Output::ExitCodes InsertEdge(const GenericEdge* edge, const unsigned int& position);
	};

	class GenericMesh : public IRotation
	{
		protected:
			vector<GenericCell*> cells; ///< Cells of the mesh
			vector<GenericFace*> faces; ///< Faces of the mesh
			vector<GenericEdge*> edges; ///< Edges of the mesh
			vector<GenericPoint*> points; ///< Points of the mesh

		public:
			GenericMesh();
			GenericMesh(const GenericMesh& mesh);
			virtual ~GenericMesh();

			const size_t NumberOfCells() const { return cells.size(); }
			const size_t NumberOfFaces() const { return faces.size(); }
			const size_t NumberOfEdges() const { return edges.size(); }
			const size_t NumberOfPoints() const { return points.size(); }
			const GenericCell* Cell(const unsigned int& position) const { return cells[position]; }
			const GenericFace* Face(const unsigned int& position) const { return faces[position]; }
			const GenericEdge* Edge(const unsigned int& position) const { return edges[position]; }
			const GenericPoint* Point(const unsigned int& position) const { return points[position]; }
			GenericCell* Cell(const unsigned int& position) { return cells[position]; }
			GenericFace* Face(const unsigned int& position) { return faces[position]; }
			GenericEdge* Edge(const unsigned int& position) { return edges[position]; }
			GenericPoint* Point(const unsigned int& position) { return points[position]; }
			const unsigned short Dimension() const { if (faces.size() > 0) return 3; else if (edges.size() > 0) return 2; else return 1; }


			/// Rotate the mesh using the rotation in IRotation
			Output::ExitCodes Rotate(const bool& inverse = false);

			virtual GenericCell* CreateCell() { return new GenericCell(cells.size()); }
			virtual GenericFace* CreateFace() { return new GenericFace(faces.size()); }
			virtual GenericEdge* CreateEdge() { return new GenericEdge(edges.size()); }
			virtual GenericPoint* CreatePoint() { return new GenericPoint(points.size()); }

			void InitializeCells(const size_t numberOfCells) { cells.reserve(numberOfCells); }
			void InitializeFaces(const size_t numberOfFaces) { faces.reserve(numberOfFaces); }
			void InitializeEdges(const size_t numberOfEdges) { edges.reserve(numberOfEdges); }
			void InitializePoints(const size_t numberOfPoints) { points.reserve(numberOfPoints); }

			Output::ExitCodes AddCell(GenericCell* cell);
			Output::ExitCodes AddFace(GenericFace* face);
			Output::ExitCodes AddEdge(GenericEdge* edge);
			Output::ExitCodes AddPoint(GenericPoint* point);

			const bool FindCoordinates(const Vector3d& coordinates, unsigned int& idPoint, const double& toll = 1.0E-5);

			const Output::ExitCodes CheckDoublePoints(const double& toll = 1.0E-7);
			const Output::ExitCodes CheckPointsInCells();
			const Output::ExitCodes CheckDoubleCells();
			const Output::ExitCodes CheckDoubleFaces();
			const Output::ExitCodes CheckDoubleEdges();
			const Output::ExitCodes CheckNeigs();

			void MovePointMesh(const unsigned int& position,const Vector3d& newCoordinate) { points[position]->SetCoordinates(newCoordinate); }
			const Output::ExitCodes CutEdgeWithPoints(const unsigned int& idEdge, const vector<Vector3d >& coordinatesPoints);

			const Output::ExitCodes UpdateFace(const unsigned int& idFace, const int& idEdge = -1);
			const Output::ExitCodes CreateFaceChild(GenericFace& face, GenericFace& faceFather, const list<unsigned int>& idEdgesFace, const list<unsigned int>& idPointFace, const bool& property = true);

			const Output::ExitCodes UpdateCell(const unsigned int& idCell, const int& idEdge = -1);
			const Output::ExitCodes CreateCellChild2D(GenericCell& cell, GenericCell& cellFather, const list<unsigned int>& idEdgesCell, const list<unsigned int>& idPointCell, const bool& property = true);

			const Output::ExitCodes ActivateFatherNodes();
			const Output::ExitCodes ActivateChildrenNodes();
			const Output::ExitCodes CleanInactiveTreeNode();
	};

	class GenericMeshImportInterface
	{
		protected:
			double maximumCellSize; ///< Size of the minimum cell of the mesh
			unsigned int minimumNumberOfCells; ///< Minimum number of cell of the mesh

			vector<int> vertexMarkers; ///< Vector of boundary conditions of domain vertices
			vector<int> edgeMarkers; ///< Vector of boundary conditions of domain edges
			vector<int> faceMarkers; ///< Vector of boundary conditions of domain faces

		public:
			GenericMeshImportInterface();
			virtual ~GenericMeshImportInterface();

			const double& MaximumCellSize() const { return maximumCellSize; }
			const unsigned int& MinimumNumberOfCells() const { return minimumNumberOfCells; }

			const vector<int>& VertexMarkers() const { return vertexMarkers; }
			const vector<int>& EdgeMarkers() const { return edgeMarkers; }
			const vector<int>& FaceMarkers() const { return faceMarkers; }

			void SetMaximumCellSize(const double& _maximumCellSize) { maximumCellSize = _maximumCellSize; }
			void SetMinimumNumberOfCells(const unsigned int& _minimumNumberOfCells) { minimumNumberOfCells = _minimumNumberOfCells; }
			void SetBoundaryConditions(const vector<int>& _vertexMarkers, const vector<int>& _edgeMarkers = vector<int>(), const vector<int>& _faceMarkers = vector<int>());

			virtual Output::ExitCodes CreateMesh(const GenericDomain& domain, GenericMesh& mesh) const = 0;
	};
}

#endif // GENERICMESH_H
