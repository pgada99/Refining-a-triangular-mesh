#include "GenericMesh.hpp"
#include "Output.hpp"
#include <set>
#define DEBUG


using namespace MainApplication;

namespace GeDiM
{
	// ***************************************************************************
	GenericTreeNode::GenericTreeNode(const unsigned int& _id)
	{
		father = NULL;

		id = _id;
		globalId = _id;
		isActive = true;

		marker = 0;
	}

	GenericTreeNode::GenericTreeNode(const GenericTreeNode& treeNode)
	{
		father = NULL;
		childs.resize(treeNode.childs.size(), NULL); ///< Child nodes in the tree.

		id = treeNode.id; ///< Internal id of the node
		globalId = treeNode.globalId; ///< Global id of the node, generally equal to id
		isActive = treeNode.isActive; ///< Tells if the node is active for system matrix computation

		marker = treeNode.marker; ///< Marker for Dirichlet and Neumann conditions
		properties = *new map<string, void*>(treeNode.properties.begin(),treeNode.properties.end());
	}
	GenericTreeNode::~GenericTreeNode()
	{
		properties.clear();
		childs.clear();
	}
	// ***************************************************************************
	Output::ExitCodes GenericTreeNode::InheritPropertiesByFather()
	{
		if(father == NULL)
			return Output::GenericError;

		const map< string, void* >& propertiesFather = father->GetAllProperties();
		if(propertiesFather.size() == 0)
			return Output::Success;

		for(map< string, void* >::const_iterator iteratorProperties = propertiesFather.begin(); iteratorProperties != propertiesFather.end(); iteratorProperties++)
			properties[iteratorProperties->first] = iteratorProperties->second;

		return Output::Success;
	}

	// ***************************************************************************
	const Output::ExitCodes GenericTreeNode::InitializeProperty(const string& key)
	{
		map<string,void*>::iterator finder;
		finder = properties.find(key);
		if(finder == properties.end())
		{
			pair <string, void*> property(key, NULL);
			properties.insert(property);
			return Output::Success;
		}
		else
			return Output::GenericError;
	}

	// ***************************************************************************
	GenericCell::GenericCell(const unsigned int& _id) : GenericTreeNode(_id)
	{
	}

	GenericCell::GenericCell(const GenericCell& cell) : GenericTreeNode(cell)
	{
		edges.resize(cell.edges.size(), NULL);
		faces.resize(cell.faces.size(), NULL);
		cells.resize(cell.cells.size(), NULL);
		points.resize(cell.points.size(), NULL);
		childs.resize(cell.childs.size(), NULL);
		father = NULL;
	}
	GenericCell::~GenericCell()
	{
		points.clear();
		edges.clear();
		faces.clear();
		cells.clear();
	}
	// ***************************************************************************
	Output::ExitCodes GenericCell::AddCell(const GenericCell* cell)
	{
		if (cell == NULL)
			return Output::GenericError;

		cells.push_back(cell);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericCell::AddFace(const GenericFace* face)
	{
		if (face == NULL)
			return Output::GenericError;

		faces.push_back(face);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericCell::AddEdge(const GenericEdge* edge)
	{
		if (edge == NULL)
			return Output::GenericError;

		edges.push_back(edge);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericCell::AddPoint(const GenericPoint* point)
	{
		if (point == NULL)
			return Output::GenericError;

		points.push_back(point);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericCell::InsertCell(const GenericCell* cell, const unsigned int& position)
	{
		if (cell == NULL || position >= cells.size())
			return Output::GenericError;

		cells[position] = cell;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericCell::InsertFace(const GenericFace* face, const unsigned int& position)
	{
		if (face == NULL || position >= faces.size())
			return Output::GenericError;

		faces[position] = face;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericCell::InsertEdge(const GenericEdge* edge, const unsigned int& position)
	{
		if (edge == NULL || position >= edges.size())
			return Output::GenericError;

		edges[position] = edge;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericCell::InsertPoint(const GenericPoint* point, const unsigned int& position)
	{
		if (point == NULL || position >= points.size())
			return Output::GenericError;

		points[position] = point;

		return Output::Success;
	}
	// ***************************************************************************
	bool GenericCell::PointInCellAndIdEdgeBoundary (const Vector3d& point, int& numEdges, const double& toll) const
	{
		for(unsigned int pntCell = 0; pntCell < NumberOfPoints(); pntCell++)
		{
			const GenericPoint* pointFirst = Point(pntCell);
			const GenericPoint* pointSecond = Point((pntCell+1)%NumberOfPoints());

			Vector3d tangentVectorEdge = pointSecond->Coordinates() - pointFirst->Coordinates();
			Vector3d tangentVectorDifference = point - pointFirst->Coordinates();

			double crossProduct = tangentVectorEdge.x() * tangentVectorDifference.y() - tangentVectorDifference.x() * tangentVectorEdge.y();


			if( crossProduct > -toll && crossProduct < toll)
			{
				if((tangentVectorEdge.x() * tangentVectorDifference.x() < -toll) || (tangentVectorEdge.y() * tangentVectorDifference.y() < -toll))
					continue;
				if(tangentVectorEdge.squaredNorm() < tangentVectorDifference.squaredNorm())
					continue;
				numEdges = pntCell;
			}

			if( crossProduct < -toll)
				return false;
		}
		return true;
	}


	// ***************************************************************************
	bool GenericCell::PointInCell (const Vector3d& point, const double& toll) const
	{
		for(unsigned int pntCell = 0; pntCell < NumberOfPoints(); pntCell++)
		{
			const GenericPoint* pointFirst = Point(pntCell);
			const GenericPoint* pointSecond = Point((pntCell+1)%NumberOfPoints());

			Vector3d tangentVectorEdge = pointSecond->Coordinates() - pointFirst->Coordinates();
			Vector3d tangentVectorDifference = point - pointFirst->Coordinates();

			double crossProduct = tangentVectorEdge.x() * tangentVectorDifference.y() - tangentVectorDifference.x() * tangentVectorEdge.y();

			if( crossProduct < -toll)
				return false;
		}
		return true;
	}

	//***************************************************************************
	// ***************************************************************************
	GenericFace::GenericFace(const unsigned int& _id) : GenericTreeNode(_id)
	{
	}

	GenericFace::GenericFace(const GenericFace& face): GenericTreeNode(face)
	{
		edges.resize(face.edges.size(), NULL);
		faces.resize(face.faces.size(), NULL);
		cells.resize(face.cells.size(), NULL);
		points.resize(face.points.size(), NULL);
		childs.resize(face.childs.size(), NULL);
		father = NULL;
	}
	GenericFace::~GenericFace()
	{
		points.clear();
		edges.clear();
		faces.clear();
		cells.clear();
	}
	// ***************************************************************************
	Output::ExitCodes GenericFace::AddCell(const GenericCell* cell)
	{
		if (cell == NULL)
			return Output::GenericError;

		cells.push_back(cell);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericFace::AddFace(const GenericFace* face)
	{
		if (face == NULL)
			return Output::GenericError;

		faces.push_back(face);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericFace::AddEdge(const GenericEdge* edge)
	{
		if (edge == NULL)
			return Output::GenericError;

		edges.push_back(edge);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericFace::AddPoint(const GenericPoint* point)
	{
		if (point == NULL)
			return Output::GenericError;

		points.push_back(point);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericFace::InsertCell(const GenericCell* cell, const unsigned int& position)
	{
		if (cell == NULL || position >= cells.size())
			return Output::GenericError;

		cells[position] = cell;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericFace::InsertFace(const GenericFace* face, const unsigned int& position)
	{
		if (face == NULL || position >= faces.size())
			return Output::GenericError;

		faces[position] = face;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericFace::InsertEdge(const GenericEdge* edge, const unsigned int& position)
	{
		if (edge == NULL || position >= edges.size())
			return Output::GenericError;

		edges[position] = edge;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericFace::InsertPoint(const GenericPoint* point, const unsigned int& position)
	{
		if (point == NULL || position >= points.size())
			return Output::GenericError;

		points[position] = point;

		return Output::Success;
	}
	// ***************************************************************************
	const Matrix3d GenericFace::RotationMatrix(double rotationTolerance)
	{
		unsigned int totalNumberVertices = points.size();
		MatrixXd Z(3, totalNumberVertices);
		MatrixXd W(3, totalNumberVertices);
		Matrix3d H;
		Vector3d V1mV0 = points[1]->Coordinates() - points[0]->Coordinates();
		Vector3d planeNormal(0.0,0.0,0.0);
		for (unsigned int i = 0; i < totalNumberVertices - 1; i++)
		{
			Vector3d VimV0 = points[i+1]->Coordinates() - points[0]->Coordinates();
			Z.col(i) = VimV0;

			double normVector = VimV0.norm();
			double angleBetweenVectors = ((VimV0 - V1mV0).norm() < rotationTolerance) ? 0.0 : acos(VimV0.dot(V1mV0) / (normVector * V1mV0.norm()));
			W.col(i) << normVector * cos(angleBetweenVectors), normVector * sin(angleBetweenVectors), 0;
			if(i > 1)
				planeNormal += ((Vector3d)(Z.col(i))).cross((Vector3d)(Z.col(0)));
		}
		planeNormal.normalize();
		Z.col(totalNumberVertices - 1) = planeNormal;
		W.col(totalNumberVertices - 1) << 0, 0, 1;
		H = W * Z.transpose();
		JacobiSVD<MatrixXd> svd(H, ComputeThinU | ComputeThinV);
		return svd.matrixV() * (svd.matrixU()).transpose();
	}
	// ***************************************************************************
	GenericEdge::GenericEdge(const unsigned int& _id) : GenericTreeNode(_id)
	{
		points.reserve(2);
	}
	GenericEdge::GenericEdge(const unsigned int& _id, const GenericPoint* origin, const GenericPoint* end) : GenericEdge(_id)
	{
		AddPoint(origin);
		AddPoint(end);
	}

	GenericEdge::GenericEdge(const GenericEdge& edge) : GenericTreeNode(edge)
	{
		edges.resize(edge.edges.size(), NULL);
		faces.resize(edge.faces.size(), NULL);
		cells.resize(edge.cells.size(), NULL);
		points.resize(2, NULL);
		childs.resize(edge.childs.size(), NULL);
		father = NULL;
	}

	GenericEdge::~GenericEdge()
	{
		points.clear();
		edges.clear();
		faces.clear();
		cells.clear();
	}
	// ***************************************************************************
	Output::ExitCodes GenericEdge::AddCell(const GenericCell* cell)
	{
		if (cell == NULL)
			return Output::GenericError;

		cells.push_back(cell);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericEdge::AddFace(const GenericFace* face)
	{
		if (face == NULL)
			return Output::GenericError;

		faces.push_back(face);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericEdge::AddEdge(const GenericEdge* edge)
	{
		if (edge == NULL)
			return Output::GenericError;

		edges.push_back(edge);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericEdge::AddPoint(const GenericPoint* point)
	{
		if (point == NULL)
			return Output::GenericError;

		points.push_back(point);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericEdge::InsertCell(const GenericCell* cell, const unsigned int& position)
	{
		if (cell == NULL || position >= cells.size())
			return Output::GenericError;

		cells[position] = cell;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericEdge::InsertFace(const GenericFace* face, const unsigned int& position)
	{
		if (face == NULL || position >= faces.size())
			return Output::GenericError;

		faces[position] = face;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericEdge::InsertEdge(const GenericEdge* edge, const unsigned int& position)
	{
		if (edge == NULL || position >= edges.size())
			return Output::GenericError;

		edges[position] = edge;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericEdge::InsertPoint(const GenericPoint* point, const unsigned int& position)
	{
		if (point == NULL || position >= points.size())
			return Output::GenericError;

		points[position] = point;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericEdge::ChangeOrientation()
	{
		const GenericPoint* tempPoint = points[0];
		if(tempPoint == NULL)
			return Output::GenericError;
		points[0] = points[1];
		points[1] = tempPoint;
		return Output::Success;
	}
	// ***************************************************************************
	GenericEdge::PositionPoint GenericEdge::PointOnEdge(const Vector3d& point, const double& toll) const
	{
		Vector3d tangentVectorEdge = Point(1)->Coordinates() - Point(0)->Coordinates();
		Vector3d tangentVectorDiff = point - Point(0)->Coordinates();

		double crossProd = tangentVectorEdge.x() * tangentVectorDiff.y() - tangentVectorDiff.x() * tangentVectorEdge.y();

		if(crossProd > toll)
			return PositionPoint::AtTheLeft;
		if(crossProd < -toll)
			return PositionPoint::AtTheRight;
		if((tangentVectorEdge.x() * tangentVectorDiff.x() < -toll) || (tangentVectorEdge.y() * tangentVectorDiff.y() < -toll))
			return PositionPoint::Behind;
		if(tangentVectorEdge.squaredNorm() < tangentVectorDiff.squaredNorm())
			return PositionPoint::Beyond;
		if(tangentVectorDiff.squaredNorm() < toll * toll)
			return PositionPoint::AtTheOrigin;
		if((point - Point(1)->Coordinates()).squaredNorm() < toll * toll)
			return PositionPoint::AtTheEnd;
		return PositionPoint::Between;
	}
	// ***************************************************************************
	GenericPoint::GenericPoint(const unsigned int& _id) : GenericTreeNode(_id)
	{
		coordinates.setZero();
	}

	GenericPoint::GenericPoint(const GenericPoint& point) : GenericTreeNode(point)
	{
		coordinates = point.coordinates;
		edges.resize(point.edges.size(), NULL);
		faces.resize(point.faces.size(), NULL);
		cells.resize(point.cells.size(), NULL);
	}

	GenericPoint::~GenericPoint()
	{
		edges.clear();
		faces.clear();
		cells.clear();
	}
	// ***************************************************************************
	Output::ExitCodes GenericPoint::SetCoordinates(const Vector3d& _coordinates)
	{
		coordinates = _coordinates;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericPoint::SetCoordinates(const double& x, const double& y, const double& z)
	{
		coordinates<< x, y, z;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericPoint::AddCell(const GenericCell* cell)
	{
		if (cell == NULL)
			return Output::GenericError;

		cells.push_back(cell);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericPoint::AddFace(const GenericFace* face)
	{
		if (face == NULL)
			return Output::GenericError;

		faces.push_back(face);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericPoint::AddEdge(const GenericEdge* edge)
	{
		if (edge == NULL)
			return Output::GenericError;

		edges.push_back(edge);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericPoint::InsertCell(const GenericCell* cell, const unsigned int& position)
	{
		if (cell == NULL || position >= cells.size())
			return Output::GenericError;

		cells[position] = cell;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericPoint::InsertFace(const GenericFace* face, const unsigned int& position)
	{
		if (face == NULL || position >= faces.size())
			return Output::GenericError;

		faces[position] = face;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericPoint::InsertEdge(const GenericEdge* edge, const unsigned int& position)
	{
		if (edge == NULL || position >= edges.size())
			return Output::GenericError;

		edges[position] = edge;

		return Output::Success;
	}


	// ***************************************************************************
	GenericMesh::GenericMesh()
	{
	}

	GenericMesh::GenericMesh(const GenericMesh& mesh) : IRotation(mesh)
	{
		unsigned int numberOfCells = mesh.cells.size();
		unsigned int numberOfFaces = mesh.faces.size();
		unsigned int numberOfEdges = mesh.edges.size();
		unsigned int numberOfPoints = mesh.points.size();

		cells.resize(numberOfCells);
		faces.resize(numberOfFaces);
		edges.resize(numberOfEdges);
		points.resize(numberOfPoints);

		for(unsigned int cel = 0; cel < numberOfCells; cel++)
			cells[cel] = new GenericCell(*mesh.cells[cel]);
		for(unsigned int fac = 0; fac < numberOfFaces; fac++)
			faces[fac] = new GenericFace(*mesh.faces[fac]);
		for(unsigned int edg = 0; edg < numberOfEdges; edg++)
			edges[edg] = new GenericEdge(*mesh.edges[edg]);
		for(unsigned int pnt = 0; pnt < numberOfPoints; pnt++)
			points[pnt] = new GenericPoint(*mesh.points[pnt]);

		for(unsigned int cel = 0; cel < numberOfCells; cel++)
		{
			GenericCell& cell = *cells[cel];
			GenericCell& cellMesh = *mesh.cells[cel];
			for(unsigned int pnt = 0; pnt < cellMesh.NumberOfPoints(); pnt++)
				cell.InsertPoint(points[cellMesh.Point(pnt)->Id()], pnt);
			for(unsigned int edg = 0; edg < cellMesh.NumberOfEdges(); edg++)
				cell.InsertEdge(edges[cellMesh.Edge(edg)->Id()], edg);
			for(unsigned int fac = 0; fac < cellMesh.NumberOfFaces(); fac++)
				cell.InsertFace(faces[cellMesh.Face(fac)->Id()], fac);
			for(unsigned int cel = 0; cel < cellMesh.NumberOfCells(); cel++)
			{
				if(cellMesh.Cell(cel) != NULL)
					cell.InsertCell(cells[cellMesh.Cell(cel)->Id()], cel);
			}

			if(cellMesh.HasFather())
				cell.SetFather(cells[cellMesh.Father()->Id()]);

			if(cellMesh.HasChilds())
			{
				for(unsigned int chd = 0; chd < cellMesh.NumberOfChilds(); chd++)
					cell.InsertChild(cells[cellMesh.Child(chd)->Id()], chd);
			}
		}

		for(unsigned int fac = 0; fac < numberOfFaces; fac++)
		{
			GenericFace& face = *faces[fac];
			GenericFace& faceMesh = *mesh.faces[fac];
			for(unsigned int pnt = 0; pnt < faceMesh.NumberOfPoints(); pnt++)
				face.InsertPoint(points[faceMesh.Point(pnt)->Id()], pnt);
			for(unsigned int edg = 0; edg < faceMesh.NumberOfEdges(); edg++)
				face.InsertEdge(edges[faceMesh.Edge(edg)->Id()], edg);
			for(unsigned int fac = 0; fac < faceMesh.NumberOfFaces(); fac++)
				face.InsertFace(faces[faceMesh.Face(fac)->Id()], fac);
			for(unsigned int cel = 0; cel < faceMesh.NumberOfCells(); cel++)
			{
				if(faceMesh.Cell(cel) != NULL)
					face.InsertCell(cells[faceMesh.Cell(cel)->Id()], cel);
			}

			if(faceMesh.HasFather())
				face.SetFather(cells[faceMesh.Father()->Id()]);

			if(faceMesh.HasChilds())
			{
				for(unsigned int chd = 0; chd < faceMesh.NumberOfChilds(); chd++)
					face.InsertChild(cells[faceMesh.Child(chd)->Id()], chd);
			}
		}

		for(unsigned int edg = 0; edg < numberOfEdges; edg++)
		{
			GenericEdge& edge = *edges[edg];
			GenericEdge& edgeMesh = *mesh.edges[edg];
			for(unsigned int pnt = 0; pnt < 2; pnt++)
				edge.InsertPoint(points[edgeMesh.Point(pnt)->Id()], pnt);
			for(unsigned int edg = 0; edg < edgeMesh.NumberOfEdges(); edg++)
				edge.InsertEdge(edges[edgeMesh.Edge(edg)->Id()], edg);
			for(unsigned int fac = 0; fac < edgeMesh.NumberOfFaces(); fac++)
				edge.InsertFace(faces[edgeMesh.Face(fac)->Id()], fac);
			for(unsigned int cel = 0; cel < edgeMesh.NumberOfCells(); cel++)
			{
				if(edgeMesh.Cell(cel) != NULL)
					edge.InsertCell(cells[edgeMesh.Cell(cel)->Id()], cel);
			}

			if(edgeMesh.HasFather())
				edge.SetFather(cells[edgeMesh.Father()->Id()]);

			if(edgeMesh.HasChilds())
			{
				for(unsigned int chd = 0; chd < edgeMesh.NumberOfChilds(); chd++)
					edge.InsertChild(cells[edgeMesh.Child(chd)->Id()], chd);
			}
		}

		for(unsigned int pnt = 0; pnt < numberOfPoints; pnt++)
		{
			GenericPoint& point = *points[pnt];
			GenericPoint& pointMesh = *mesh.points[pnt];
			for(unsigned int edg = 0; edg < pointMesh.NumberOfEdges(); edg++)
				point.InsertEdge(edges[pointMesh.Edge(edg)->Id()], edg);
			for(unsigned int fac = 0; fac < pointMesh.NumberOfFaces(); fac++)
				point.InsertFace(faces[pointMesh.Face(fac)->Id()], fac);
			for(unsigned int cel = 0; cel < pointMesh.NumberOfCells(); cel++)
				point.InsertCell(cells[pointMesh.Cell(cel)->Id()], cel);
		}
	}

	GenericMesh::~GenericMesh()
	{
		for (vector<GenericPoint*>::iterator pointPtr = points.begin(); pointPtr != points.end(); pointPtr++)
			delete *pointPtr;

		for (vector<GenericEdge*>::iterator edgePtr = edges.begin(); edgePtr != edges.end(); edgePtr++)
			delete *edgePtr;

		for (vector<GenericFace*>::iterator facePtr = faces.begin(); facePtr != faces.end(); facePtr++)
			delete *facePtr;

		for (vector<GenericCell*>::iterator cellPtr = cells.begin(); cellPtr != cells.end(); cellPtr++)
			delete *cellPtr;

		points.clear();
		edges.clear();
		faces.clear();
		cells.clear();
	}
	// ***************************************************************************
	Output::ExitCodes GenericMesh::Rotate(const bool& inverse)
	{
		for (vector<GenericPoint*>::iterator pointPtr = points.begin(); pointPtr != points.end(); pointPtr++)
		{
			GenericPoint& point = **pointPtr;
			Vector3d& coordinates = point.Coordinates();

			coordinates = RotatePoint(coordinates, true, inverse);
		}

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericMesh::AddCell(GenericCell* cell)
	{
		if (cell == NULL)
			return Output::GenericError;

		cells.push_back(cell);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericMesh::AddFace(GenericFace* face)
	{
		if (face == NULL)
			return Output::GenericError;

		faces.push_back(face);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericMesh::AddEdge(GenericEdge* edge)
	{
		if (edge == NULL)
			return Output::GenericError;

		edges.push_back(edge);

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes GenericMesh::AddPoint(GenericPoint* point)
	{
		if (point == NULL)
			return Output::GenericError;

		points.push_back(point);

		return Output::Success;
	}
	// ***************************************************************************
	const bool GenericMesh::FindCoordinates(const Vector3d& coordinates, unsigned int& idPoint, const double& toll)
	{
		double squaredToll = toll* toll;
		for(unsigned int numPnt = 0; numPnt < points.size(); numPnt++)
		{
			const GenericPoint& point = *points[numPnt];
			Vector3d diff = coordinates - point.Coordinates();
			if(diff.squaredNorm() < squaredToll)
			{
				idPoint = point.Id();
				return true;
			}
		}
		return false;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericMesh::CheckDoublePoints(const double& toll)
	{
		double squaredToll = toll*toll;
		for(unsigned int pnt = 0; pnt < points.size() - 1; pnt++)
		{
			GenericPoint& point = *points[pnt];
			for(unsigned int pnt2 = pnt+1; pnt2 < points.size(); pnt2++)
			{
				GenericPoint& point2 = *points[pnt2];
				Vector3d diff = point2.Coordinates() - point.Coordinates();
				if(diff.squaredNorm() < squaredToll)
				{
					Output::PrintWarningMessage("Point %d and Point %d have the same coordinates", true, point.Id(), point2.Id());

				}
			}
		}
		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericMesh::CheckPointsInCells()
	{
		set<unsigned int> pointsInCells;
		for(unsigned int cel = 0; cel < cells.size(); cel++)
		{
			GenericCell& cell = *cells[cel];
			for(unsigned int pnt = 0; pnt < cell.NumberOfPoints(); pnt++)
				pointsInCells.insert(cell.Point(pnt)->Id());
		}
		if(pointsInCells.size() != points.size())
		{
			Output::PrintErrorMessage("There are %d point/s not used ", true, pointsInCells.size() - points.size());
			return Output::GenericError;
		}

		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericMesh::CheckDoubleCells()
	{
		set<unsigned int> firstCell;
		for(unsigned int cel = 0; cel < cells.size()-1; cel++)
		{
			GenericCell& cell = *cells[cel];
			unsigned int numFacesOrEdges = cell.NumberOfFaces();
			bool controlEdges = false;
			if(numFacesOrEdges == 0)
			{
				controlEdges = true;
				numFacesOrEdges = cell.NumberOfEdges();
			}
			if(!controlEdges)
			{
				if(numFacesOrEdges < 4)
				{
					Output::PrintErrorMessage("In cell %d there are %d faces ", true, cell.Id(), numFacesOrEdges);
					return Output::GenericError;
				}
			}
			else
			{
				if(numFacesOrEdges < 3)
				{
					Output::PrintErrorMessage("In cell %d there are %d edges ", true, cell.Id(), numFacesOrEdges);
					return Output::GenericError;
				}
			}

			if(!controlEdges)
				for(unsigned int fac = 0; fac < numFacesOrEdges; fac++)
					firstCell.insert(cell.Face(fac)->Id());
			else
				for(unsigned int edg = 0; edg < numFacesOrEdges; edg++)
					firstCell.insert(cell.Edge(edg)->Id());

			if(firstCell.size() != numFacesOrEdges)
			{
				if(!controlEdges)
					Output::PrintErrorMessage("In cell %d there are doubled faces", true, cell.Id());
				else
					Output::PrintErrorMessage("In cell %d there are doubled edges", true, cell.Id());
				return Output::GenericError;
			}

			set<unsigned int> secondCell;
			for(unsigned int cel2 = cel+1; cel2 < cells.size(); cel2++)
			{
				GenericCell& cell2 = *cells[cel2];
				if(cell2.NumberOfFaces() != numFacesOrEdges)
					continue;

				if(!controlEdges)
					for(unsigned int fac2 = 0; fac2 < numFacesOrEdges; fac2++)
						secondCell.insert(cell2.Face(fac2)->Id());
				else
					for(unsigned int edg2 = 0; edg2 < numFacesOrEdges; edg2++)
						secondCell.insert(cell2.Edge(edg2)->Id());

				if(secondCell.size() != numFacesOrEdges)
				{
					if(!controlEdges)
						Output::PrintErrorMessage("In cell %d there are doubled faces ", true, cell2.Id());
					else
						Output::PrintErrorMessage("In cell %d there are doubled edges ", true, cell2.Id());
					return Output::GenericError;
				}

				set<unsigned int>::iterator it = firstCell.begin();
				set<unsigned int>::iterator it2 = secondCell.begin();
				bool cont = true;
				for(unsigned int i = 0; i < numFacesOrEdges; i++)
				{
					if(*it != *it2)
					{
						cont = false;
						break;
					}
					it++;
					it2++;
				}
				if(cont)
				{
					Output::PrintErrorMessage("There is a doubled cell with id %d and %d ", true, cell.Id(), cell2.Id());
					return Output::GenericError;
				}
				secondCell.clear();
			}
			firstCell.clear();
		}
		return Output::Success;
	}

	// ***************************************************************************
	const Output::ExitCodes GenericMesh::CheckDoubleFaces()
	{
		set<unsigned int> firstFace;
		for(unsigned int fac = 0; fac < faces.size()-1; fac++)
		{
			GenericFace& face = *faces[fac];
			unsigned int numEdges = face.NumberOfEdges();
			if(numEdges < 3)
			{
				Output::PrintErrorMessage("In face %d there are %d edges ", true, face.Id(), numEdges);
				return Output::GenericError;
			}

			for(unsigned int edg = 0; edg < numEdges; edg++)
				firstFace.insert(face.Edge(edg)->Id());

			if(firstFace.size() != numEdges)
			{
				Output::PrintErrorMessage("In face %d there are doubled edges", true, face.Id());
				return Output::GenericError;
			}

			set<unsigned int> secondFace;
			for(unsigned int fac2 = fac+1; fac2 < faces.size(); fac2++)
			{
				GenericFace& face2 = *faces[fac2];
				if(face2.NumberOfFaces() != numEdges)
					continue;

				for(unsigned int edg2 = 0; edg2 < numEdges; edg2++)
					secondFace.insert(face2.Edge(edg2)->Id());

				if(secondFace.size() != numEdges)
				{
					Output::PrintErrorMessage("In face %d there are doubled edges ", true, face2.Id());
					return Output::GenericError;
				}

				set<unsigned int>::iterator it = firstFace.begin();
				set<unsigned int>::iterator it2 = secondFace.begin();
				bool cont = true;
				for(unsigned int i = 0; i < numEdges; i++)
				{
					if(*it != *it2)
					{
						cont = false;
						break;
					}
					it++;
					it2++;
				}
				if(cont)
				{
					Output::PrintErrorMessage("There is a doubled face with id %d and %d ", true, face.Id(), face2.Id());
					return Output::GenericError;
				}
				secondFace.clear();
			}
			firstFace.clear();
		}
		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericMesh::CheckDoubleEdges()
	{
		set<unsigned int> firstEdge;
		for(unsigned int edg = 0; edg < edges.size()-1; edg++)
		{
			GenericEdge& edge = *edges[edg];
			unsigned int numPoints = edge.NumberOfPoints();
			if(numPoints != 2)
			{
				Output::PrintErrorMessage("In edge %d there are %d points ", true, edge.Id(), numPoints);
				return Output::GenericError;
			}

			for(unsigned int edg = 0; edg < 2; edg++)
				firstEdge.insert(edge.Point(edg)->Id());

			if(firstEdge.size() != 2)
			{
				Output::PrintErrorMessage("In edge %d there are doubled points", true, edge.Id());
				return Output::GenericError;
			}

			set<unsigned int> secondEdge;
			for(unsigned int edg2 = edg+1; edg2 < edges.size(); edg2++)
			{
				GenericEdge& edge2 = *edges[edg2];

				for(unsigned int edg2 = 0; edg2 < 2; edg2++)
					secondEdge.insert(edge2.Point(edg2)->Id());

				if(secondEdge.size() != numPoints)
				{
					Output::PrintErrorMessage("In edge %d there are doubled point ", true, edge2.Id());
					return Output::GenericError;
				}

				set<unsigned int>::iterator it = firstEdge.begin();
				set<unsigned int>::iterator it2 = secondEdge.begin();
				bool cont = true;
				for(unsigned int i = 0; i < numPoints; i++)
				{
					if(*it != *it2)
					{
						cont = false;
						break;
					}
					it++;
					it2++;
				}
				if(cont)
				{
					Output::PrintErrorMessage("There is a doubled edge with id %d and %d ", true, edge.Id(), edge2.Id());
					return Output::GenericError;
				}
				secondEdge.clear();
			}
			firstEdge.clear();
		}
		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericMesh::CheckNeigs()
	{
		for(unsigned int cel = 0; cel < cells.size()-1; cel++)
		{
			GenericCell& cell = *cells[cel];
			unsigned int idCell = cell.Id();
			unsigned int idCellToControl = 0;
			unsigned int dimension = Dimension();
			switch(dimension)
			{
				case 3:
				{
					for(unsigned int fac = 0; fac < cell.NumberOfFaces(); fac++)
					{
						const GenericFace& face = *cell.Face(fac);
						if(face.Cell(0) == NULL )
						{
							if(cell.Cell(fac) == NULL)
								continue;
							else
							{
								Output::PrintErrorMessage("Position neigs not right", true);
								return Output::GenericError;
							}
						}
						else if(face.Cell(0) != NULL && face.Cell(0)->Id() != idCell)
						{
							idCellToControl = face.Cell(0)->Id();
							if(cell.Cell(fac)->Id() == idCellToControl)
								continue;
							else
							{
								Output::PrintErrorMessage("Position neigs not right", true);
								return Output::GenericError;
							}
						}
						if(face.Cell(1) == NULL )
						{
							if(cell.Cell(fac) == NULL)
								continue;
							else
							{
								Output::PrintErrorMessage("Position neigs not right", true);
								return Output::GenericError;
							}
						}
						else if(face.Cell(1) != NULL && face.Cell(1)->Id() != idCell)
						{
							idCellToControl = face.Cell(1)->Id();
							if(cell.Cell(fac)->Id() == idCellToControl)
								continue;
							else
							{
								Output::PrintErrorMessage("Position neigs not right", true);
								return Output::GenericError;
							}
						}
					}

				}
					break;
				case 2:
				{
					for(unsigned int edg = 0; edg < cell.NumberOfEdges(); edg++)
					{
						const GenericEdge& edge = *cell.Edge(edg);
						if(edge.Cell(0) == NULL )
						{
							if(cell.Cell(edg) == NULL)
								continue;
							else
							{
								Output::PrintErrorMessage("Position neigs not right", true);
								return Output::GenericError;
							}
						}
						else if(edge.Cell(0) != NULL && edge.Cell(0)->Id() != idCell)
						{
							idCellToControl = edge.Cell(0)->Id();
							if(cell.Cell(edg)->Id() == idCellToControl)
								continue;
							else
							{
								Output::PrintErrorMessage("Position neigs not right", true);
								return Output::GenericError;
							}
						}
						if(edge.Cell(1) == NULL )
						{
							if(cell.Cell(edg) == NULL)
								continue;
							else
							{
								Output::PrintErrorMessage("Position neigs not right", true);
								return Output::GenericError;
							}
						}
						else if(edge.Cell(1) != NULL && edge.Cell(1)->Id() != idCell)
						{
							idCellToControl = edge.Cell(1)->Id();
							if(cell.Cell(edg)->Id() == idCellToControl)
								continue;
							else
							{
								Output::PrintErrorMessage("Position neigs not right", true);
								return Output::GenericError;
							}
						}
					}
				}
					break;
				case 1:
				{
					for(unsigned int pnt = 0; pnt < cell.NumberOfPoints(); pnt++)
					{
						const GenericPoint& point = *cell.Point(pnt);
						if(point.Cell(0) == NULL )
						{
							if(cell.Cell(pnt) == NULL)
								continue;
							else
							{
								Output::PrintErrorMessage("Position neigs not right", true);
								return Output::GenericError;
							}
						}
						else if(point.Cell(0) != NULL && point.Cell(0)->Id() != idCell)
						{
							idCellToControl = point.Cell(0)->Id();
							if(cell.Cell(pnt)->Id() == idCellToControl)
								continue;
							else
							{
								Output::PrintErrorMessage("Position neigs not right", true);
								return Output::GenericError;
							}
						}
						if(point.Cell(1) == NULL )
						{
							if(cell.Cell(pnt) == NULL)
								continue;
							else
							{
								Output::PrintErrorMessage("Position neigs not right", true);
								return Output::GenericError;
							}
						}
						else if(point.Cell(1) != NULL && point.Cell(1)->Id() != idCell)
						{
							idCellToControl = point.Cell(1)->Id();
							if(cell.Cell(pnt)->Id() == idCellToControl)
								continue;
							else
							{
								Output::PrintErrorMessage("Position neigs not right", true);
								return Output::GenericError;
							}
						}
					}
				}
					break;
				default:
					break;
			}
		}
		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericMesh::CutEdgeWithPoints(const unsigned int& idEdge, const vector<Vector3d>& coordinatesPoints)
	{
		unsigned int numPoints = coordinatesPoints.size();
		unsigned int numEdges = numPoints +1;

		vector< GenericPoint*> pointTraces(numPoints, NULL);
		vector<GenericEdge*> newEdges(numEdges, NULL);

		GenericEdge& edge = *edges[idEdge];
		edge.SetState(false);

		GenericPoint* firstPoint = points[edge.Point(0)->Id()];
		GenericPoint* secondPoint = points[edge.Point(1)->Id()];

		unsigned int marker = edge.Marker();
		for(unsigned int pnt = 0; pnt < numPoints; pnt++)
		{
			pointTraces[pnt] = CreatePoint();
			pointTraces[pnt]->SetMarker(marker);
			AddPoint(pointTraces[pnt]);
			pointTraces[pnt]->InitializeEdges(2);
			pointTraces[pnt]->SetCoordinates(coordinatesPoints[pnt]);
		}

		unsigned int numberCell = edge.NumberOfCells();

		edge.InitializeChilds(numEdges);
		for(unsigned int edg = 0; edg < numEdges; edg++)
		{
			newEdges[edg] = CreateEdge();
			AddEdge(newEdges[edg]);
			newEdges[edg]->SetFather(&edge);
			newEdges[edg]->AllocateCells(numberCell);
			newEdges[edg]->SetMarker(marker);
			edge.AddChild(newEdges[edg]);
		}

		newEdges[0]->AddPoint(firstPoint);
		newEdges[0]->AddPoint(pointTraces[0]);
		newEdges[numPoints]->AddPoint(pointTraces[numPoints-1]);
		newEdges[numPoints]->AddPoint(secondPoint);

		for(unsigned int edg = 0; edg < numPoints - 1; edg++)
		{
			newEdges[edg+1]->AddPoint(pointTraces[edg]);
			newEdges[edg+1]->AddPoint(pointTraces[edg+1]);
		}

		for(unsigned int pnt = 0; pnt < numPoints; pnt++)
		{
			pointTraces[pnt]->AddEdge(newEdges[pnt]);
			pointTraces[pnt]->AddEdge(newEdges[pnt+1]);
		}

		for(unsigned int edg = 0; edg < firstPoint->NumberOfEdges(); edg++)
			if(firstPoint->Edge(edg)->Id() == edge.Id())
				firstPoint->InsertEdge(newEdges[0],edg);
		for(unsigned int edg = 0; edg < secondPoint->NumberOfEdges(); edg++)
			if(secondPoint->Edge(edg)->Id() == edge.Id())
				secondPoint->InsertEdge(newEdges[numEdges-1],edg);

		if(edge.Cell(0) != NULL)
		{
			for(unsigned int edg = 0; edg < numEdges; edg++)
				newEdges[edg]->InsertCell(edge.Cell(0),0);
			for(unsigned int pnt = 0; pnt < numPoints; pnt++)
				pointTraces[pnt]->AddCell(edge.Cell(0));

			UpdateCell(edge.Cell(0)->Id());
		}
		if(edge.Cell(1) != NULL)
		{
			for(unsigned int edg = 0; edg < numEdges; edg++)
				newEdges[edg]->InsertCell(edge.Cell(1),1);
			for(unsigned int pnt = 0; pnt < numPoints; pnt++)
				pointTraces[pnt]->AddCell(edge.Cell(1));

			UpdateCell(edge.Cell(1)->Id());
		}
		return Output::Success;
	}

	// ***************************************************************************
	const Output::ExitCodes GenericMesh::UpdateCell(const unsigned int& idCell, const int& idEdge)
	{
		GenericCell& cell = *cells[idCell];
		if(!cell.IsActive())
			return Output::Success;
		GenericCell& cellChild = *CreateCell();

		cell.InitializeChilds(1);
		cell.AddChild(&cellChild);
		cell.SetState(false);
		cellChild.SetFather(&cell);
		cellChild.InheritPropertiesByFather();

		unsigned int dimension = Dimension();
		switch(dimension)
		{
			case 2:
			{

				unsigned int numberEdges = cell.NumberOfEdges(); //equal number points
				unsigned int numberEdgesChild = 0;

				vector<unsigned int> pointsIdCutEdge(2,0); //id dei punti del lato tagliato
				vector<unsigned int> newIdPoint; // nuovo punto da aggiungere

				if(idEdge == -1)
				{
					for(unsigned int numEd = 0; numEd < numberEdges; numEd++)
					{
						unsigned idEdgeCut = cell.Edge(numEd)->Id();
						GenericEdge* edge = edges[idEdgeCut];
						if(edge->HasChilds())
						{
							pointsIdCutEdge[0] = edge->Point(0)->Id();
							pointsIdCutEdge[1] = edge->Point(1)->Id();
							numberEdgesChild = edge->NumberOfChilds();
							break;
						}
					}
				}
				else
				{
					GenericEdge* edge = edges[idEdge];
					pointsIdCutEdge[0] = edge->Point(0)->Id();
					pointsIdCutEdge[1] = edge->Point(1)->Id();
					numberEdgesChild = edge->NumberOfChilds();
				}

				for(unsigned int numPoint = 0; numPoint < numberEdges; numPoint++)
				{
					unsigned int idPoint = cell.Point(numPoint)->Id();
					unsigned int idPointNext = cell.Point((numPoint+1)%numberEdges)->Id();
					if(idPoint == pointsIdCutEdge[1] && idPointNext == pointsIdCutEdge[0])
					{
						pointsIdCutEdge[0] = idPoint;
						pointsIdCutEdge[1] = idPointNext;
						break;
					}
				}

				unsigned int numEdgesCellChild = numberEdges - 1 + numberEdgesChild;
				cellChild.InitializeEdges(numEdgesCellChild);
				cellChild.InitializePoints(numEdgesCellChild);
				cellChild.AllocateCells(numEdgesCellChild);
				AddCell(&cellChild);

				newIdPoint.reserve(numEdgesCellChild-1);

				//CICLO SUI LATI DELLA CELLA
				for(unsigned int numEd = 0; numEd < numberEdges; numEd++)
				{
					unsigned int idEdge = cell.Edge(numEd)->Id();
					GenericEdge* edge = edges[idEdge];
					//SE NON HA FIGLI IL LATO:
					//1) AGGIUNGO IL LATO ALLA CELLA
					//2) AGGIORNO LA CELLA DEL LATO PADRE CON LA CELLA FIGLIO
					if(!edge->HasChilds())
					{
						cellChild.AddEdge(edge);
						for(unsigned int neigCell = 0; neigCell < 2; neigCell++)
							if(edge->Cell(neigCell) != NULL)
								if(edge->Cell(neigCell)->Id() == cell.Id())
									edge->InsertCell(&cellChild, neigCell);
					}
					//SE HA FIGLI IL LATO:
					//1) SALVO GLI ID DEI PUNTI DEL LATO PADRE CHE E' STATO TAGLIATO
					//2) AGGIUNGO ALLA CELLA I FIGLI DEL LATO
					//3) AI FIGLI DEL LATO AGGIORNO LA CELLA PADRE CON LA CELLA FIGLIO
					//4) TROVO L'ID DEL NUOVO PUNTO DA AGGIUNGERE
					//5) AGGIORNO GLI ID DEI PUNTI DEL LATO PADRE CON QUELLI DEL FIGLIO TEMPORANEAMENTE
					//   PER TROVARE L'EVENTUALE NUOVO ID DEL PUNTO SE NE SONO DUE
					//6) SALVO GLI ID DEI PUNTI DEL LATO PADRE CHE E' STATO TAGLIATO
					else
					{
						bool inversePushBack = false;
						unsigned int idEdge = edge->Child(0)->Id();
						GenericEdge* edgeChild = edges[idEdge];
						if(edgeChild->Point(0)->Id() == pointsIdCutEdge[1])
							inversePushBack = true;

						unsigned int counter = 0;
						if(inversePushBack)
						{
							for(int numChild = edge->NumberOfChilds() - 1; numChild >= 0 ; numChild--)
							{
								unsigned int idEdge = edge->Child(numChild)->Id();
								GenericEdge* edgeChild = edges[idEdge];
								cellChild.AddEdge(edgeChild);

								for(unsigned int neigCell = 0; neigCell < 2; neigCell++)
									if(edge->Cell(neigCell) != NULL)
										if(edge->Cell(neigCell)->Id() == cell.Id())
											edgeChild->InsertCell(&cellChild, neigCell);

								if(counter < edge->NumberOfChilds() - 1)
								{
									if(edgeChild->Point(0)->Id() != pointsIdCutEdge[0] && edgeChild->Point(0)->Id() != pointsIdCutEdge[1])
										newIdPoint.push_back(edgeChild->Point(0)->Id());
									else
										newIdPoint.push_back(edgeChild->Point(1)->Id());
								}
								pointsIdCutEdge[0] = edgeChild->Point(0)->Id();
								pointsIdCutEdge[1] = edgeChild->Point(1)->Id();
								counter++;
							}
							pointsIdCutEdge[0] = edge->Point(0)->Id();
							pointsIdCutEdge[1] = edge->Point(1)->Id();
						}
						else
						{
							for(unsigned int numChild =  0; numChild < edge->NumberOfChilds() ; numChild++)
							{
								unsigned int idEdge = edge->Child(numChild)->Id();
								GenericEdge* edgeChild = edges[idEdge];
								cellChild.AddEdge(edgeChild);

								for(unsigned int neigCell = 0; neigCell < 2; neigCell++)
									if(edge->Cell(neigCell) != NULL)
										if(edge->Cell(neigCell)->Id() == cell.Id())
											edgeChild->InsertCell(&cellChild, neigCell);

								if(counter < edge->NumberOfChilds() - 1)
								{
									if(edgeChild->Point(0)->Id() != pointsIdCutEdge[0] && edgeChild->Point(0)->Id() != pointsIdCutEdge[1])
										newIdPoint.push_back(edgeChild->Point(0)->Id());
									else
										newIdPoint.push_back(edgeChild->Point(1)->Id());
								}
								pointsIdCutEdge[0] = edgeChild->Point(0)->Id();
								pointsIdCutEdge[1] = edgeChild->Point(1)->Id();
								counter++;
							}
						}
						pointsIdCutEdge[0] = edge->Point(1)->Id();
						pointsIdCutEdge[1] = edge->Point(0)->Id();
					}
				}


				for(unsigned int numPoint = 0; numPoint < numberEdges; numPoint++)
				{
					unsigned int idPoint = cell.Point(numPoint)->Id();
					unsigned int idPointNext = cell.Point((numPoint+1)%numberEdges)->Id();
					GenericPoint* point = points[idPoint];
					cellChild.AddPoint(point);
					for(unsigned int cellPosition = 0; cellPosition < point->NumberOfCells(); cellPosition++)
					{
						//Tolgo il padre dai punti ed inserisco il figlio
						if(cell.Id() == point->Cell(cellPosition)->Id())
							point->InsertCell(&cellChild,cellPosition);
					}
					if((idPoint == pointsIdCutEdge[0] && idPointNext == pointsIdCutEdge[1]) || (idPoint == pointsIdCutEdge[1] && idPointNext == pointsIdCutEdge[0]))
					{
						for(unsigned int newId = 0; newId < newIdPoint.size() ; newId++)
						{
							GenericPoint* newPoint = points[newIdPoint[newId]];
							cellChild.AddPoint(newPoint);
							for(unsigned int cellPosition = 0; cellPosition < newPoint->NumberOfCells(); cellPosition++)
							{
								//Tolgo il padre dai punti ed inserisco il figlio
								if(cell.Id() == newPoint->Cell(cellPosition)->Id())
									newPoint->InsertCell(&cellChild,cellPosition);
							}
						}
					}
				}

				for(unsigned int numEd = 0; numEd < numEdgesCellChild; numEd++)
				{
					const GenericEdge& edge = *cellChild.Edge(numEd);
					for(unsigned int numNeig = 0; numNeig < edge.NumberOfCells(); numNeig++)
					{
						if(edge.Cell(numNeig) != NULL)
						{
							GenericCell& cellNeigh = *cells[edge.Cell(numNeig)->Id()];
							if(cellNeigh.Id() != cellChild.Id())
								cellChild.InsertCell(&cellNeigh, numEd);
							for(unsigned int neig = 0; neig < cellNeigh.NumberOfCells(); neig++)
								if(cellNeigh.Cell(neig) != NULL && cell.Id() == cellNeigh.Cell(neig)->Id())
									cellNeigh.InsertCell(&cellChild, neig);
						}
					}
				}
			}
				break;
			case 3:
			{
				set<unsigned int> pointsSetChild;
				set<unsigned int> edgesSetChild;
				set<unsigned int> facesSetChild;

				//Adding points to Cell Child
				for(unsigned int pt = 0; pt < cell.NumberOfPoints(); pt++)
					pointsSetChild.insert(cell.Point(pt)->Id());

				//Adding edges to Cell Child
				for(unsigned int edg = 0; edg < cell.NumberOfEdges(); edg++)
				{
					const GenericEdge& edge = *cell.Edge(edg);
					if(!edge.HasChilds())
					{
						edgesSetChild.insert(edge.Id());
					}
					else
					{
						for(unsigned int edgeChild = 0; edgeChild < edge.NumberOfChilds(); edgeChild++)
							edgesSetChild.insert(edge.Child(edgeChild)->Id());
					}
				}

				for(unsigned int fac = 0; fac < cell.NumberOfFaces(); fac++)
				{
					//in caso di cvx devono essere le ultime foglie dell'albero
					const GenericFace& face = *cell.Face(fac);
					if(!face.HasChilds())
						facesSetChild.insert(face.Id());
					else
					{
						for(unsigned int faceChild = 0; faceChild < face.NumberOfChilds(); faceChild++)
							facesSetChild.insert(face.Child(faceChild)->Id());
						GenericFace& faceChildReference = *faces[face.Child(0)->Id()];
						for(unsigned int pt = 0; pt < faceChildReference.NumberOfPoints(); pt++)
							pointsSetChild.insert(faceChildReference.Point(pt)->Id());

						for(unsigned int pt = 0; pt < faceChildReference.NumberOfEdges(); pt++)
							edgesSetChild.insert(faceChildReference.Edge(pt)->Id());

					}
				}


				cellChild.AllocateCells(facesSetChild.size());
				cellChild.InitializeFaces(facesSetChild.size());
				cellChild.InitializeEdges(edgesSetChild.size());
				cellChild.InitializePoints(pointsSetChild.size());
				AddCell(&cellChild);

				for(set<unsigned int>::iterator it=pointsSetChild.begin(); it!=pointsSetChild.end(); ++it)
				{
					cellChild.AddPoint(points[*it]);
					for(unsigned int pointCell = 0; pointCell < points[*it]->NumberOfCells(); pointCell++)
					{
						if(points[*it]->Cell(pointCell)->Id() == cell.Id())
							points[*it]->InsertCell(&cellChild, pointCell);
					}
				}

				for(set<unsigned int>::iterator it = edgesSetChild.begin(); it != edgesSetChild.end(); ++it)
				{
					cellChild.AddEdge(edges[*it]);
					if(!edges[*it]->HasChilds())
					{
						for(unsigned int edgeCell = 0; edgeCell < edges[*it]->NumberOfCells(); edgeCell++)
						{
							if(edges[*it]->Cell(edgeCell)->Id() == cell.Id())
								edges[*it]->InsertCell(&cellChild,edgeCell);
						}
					}
					else
					{
						edges[*it]->AddCell(&cellChild);
					}
				}


				unsigned int positionNeigCell = 0;
				for(set<unsigned int>::iterator fc = facesSetChild.begin(); fc != facesSetChild.end(); ++fc)
				{
					GenericFace& faceReference = *faces[*fc];

					if(!faceReference.HasChilds())
					{
						cellChild.AddFace(&faceReference);
						for(unsigned int neigCell = 0; neigCell < 2; neigCell++)
						{
							if(faceReference.Cell(neigCell) != NULL)
							{
								GenericCell& cellNeigh = *cells[faceReference.Cell(neigCell)->Id()];
								if(cellNeigh.Id() == cell.Id())
									faceReference.InsertCell(&cellChild, neigCell);
								else
									cellChild.InsertCell(&cellNeigh, positionNeigCell);

								for(unsigned int neig = 0; neig < cellNeigh.NumberOfCells(); neig++)
									if(cellNeigh.Cell(neig) != NULL && cell.Id() == cellNeigh.Cell(neig)->Id())
										cellNeigh.InsertCell(&cellChild, neig);
							}
						}
						positionNeigCell++;
					}
					else
					{
						for(unsigned int fcChild = 0; fcChild < faceReference.NumberOfChilds(); fcChild++)
						{
							GenericFace& faceChildReference = *faces[faceReference.Child(fcChild)->Id()];
							cellChild.AddFace(&faceChildReference);
							//updating cells of the child faces of the face in input
							for(unsigned int neigCell = 0; neigCell < 2; neigCell++)
							{
								if(faceReference.Cell(neigCell) != NULL)
								{
									GenericCell& cellNeigh = *cells[faceReference.Cell(neigCell)->Id()];
									if(cellNeigh.Id() == cell.Id())
										faceReference.InsertCell(&cellChild, neigCell);
									else
										cellChild.InsertCell(faceReference.Cell(neigCell), positionNeigCell);

									for(unsigned int neig = 0; neig < cellNeigh.NumberOfCells(); neig++)
										if(cellNeigh.Cell(neig) != NULL && cell.Id() == cellNeigh.Cell(neig)->Id())
											cellNeigh.InsertCell(&cellChild, neig);
								}
							}
							positionNeigCell++;
						}
					}
				}
			}
				break;
			default:
				break;
		}
		return Output::Success;
	}

	// ***************************************************************************
	const Output::ExitCodes GenericMesh::CreateCellChild2D(GenericCell& cell, GenericCell& cellFather, const list<unsigned int>& idEdgesCell, const list<unsigned int>& idPointCell, const bool& property)
	{
		//CONTROLLARE CHE AI LATI CI SIA TUTTO
		cellFather.AddChild(&cell);
		cellFather.SetState(false);

		cell.SetFather(&cellFather);
		unsigned int numEdgesCell = idEdgesCell.size();
		cell.InitializePoints(numEdgesCell);
		cell.AllocateEdges(numEdgesCell);
		cell.AllocateCells(numEdgesCell);
		if(property)
			cell.InheritPropertiesByFather();

		//CICLO SUGLI ID DEI PUNTI DA AGGIUNGERE ALLA PAGINA
		//1) Aggiungo il punto alla cella
		//2) Aggiorno la cella padre con la cella figlio ed inserisco update = true
		//3) Se Update rimane false significa che il punto è nuovo e non ha la cella padre e bisogna
		//   aggiungere direttamente il figlio
		for(list<unsigned int>::const_iterator iteratorId = idPointCell.begin(); iteratorId != idPointCell.end(); iteratorId++)
		{
			GenericPoint* point = points[*iteratorId];
			cell.AddPoint(point);
			bool update = false;
			for(unsigned int cellPosition = 0; cellPosition < point->NumberOfCells(); cellPosition++)
			{
				//Tolgo il padre dai punti ed inserisco il figlio
				if(cellFather.Id() == point->Cell(cellPosition)->Id())
				{
					update = true;
					point->InsertCell(&cell, cellPosition);
					break;
				}
			}
			if(!update)
				point->AddCell(&cell);
		}

		unsigned int numEdge = 0;
		for(list<unsigned int>::const_iterator iteratorId = idEdgesCell.begin(); iteratorId != idEdgesCell.end(); iteratorId++)
		{
			GenericEdge* edge = edges[*iteratorId];
			cell.InsertEdge(edge, numEdge);

			//CICLO SULLE CELLE VICINE AL LATO:
			//1) AGGIORNARE LA CELLA DEL LATO CON LA CELLA FIGLIO
			//2) AGGIUNGO LA CELLA VICINA ALLA CELLA CHE STO CREANDO
			for(unsigned int neigCellEdge = 0;  neigCellEdge < 2; neigCellEdge++)
			{
				if(edge->Cell(neigCellEdge) != NULL)
				{
					const unsigned int& idNeigCell = edge->Cell(neigCellEdge)->Id();
					if(idNeigCell == cellFather.Id())
						edge->InsertCell(&cell,neigCellEdge);
					else if((idNeigCell != cellFather.Id()) && (idNeigCell != cell.Id()))
					{
						cell.InsertCell(cells[idNeigCell], numEdge);
						GenericCell& cellNeigh = *cells[idNeigCell];
						for(unsigned int numNeigh = 0; numNeigh < cellNeigh.NumberOfCells(); numNeigh++)
							if(cellNeigh.Cell(numNeigh) != NULL && cellNeigh.Cell(numNeigh)->Id() == cellFather.Id())
								cellNeigh.InsertCell(&cell, numNeigh);
					}
				}
			}
			numEdge++;
		}
		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericMesh::UpdateFace(const unsigned int& idFace, const int& idEdge)
	{
		GenericFace& face = *faces[idFace];
		GenericFace& faceChild = *CreateFace();

		face.InitializeChilds(1);
		face.AddChild(&faceChild);
		face.SetState(false);
		faceChild.SetFather(&face);
		faceChild.InheritPropertiesByFather();

		unsigned int numberEdges = face.NumberOfEdges(); //equal number points
		unsigned int numberEdgesChild = 0;

		vector<unsigned int> pointsIdCutEdge(2,0); //id dei punti del lato tagliato
		vector<unsigned int> newIdPoint; // nuovo punto da aggiungere

		if(idEdge == -1)
		{
			for(unsigned int numEd = 0; numEd < numberEdges; numEd++)
			{
				unsigned idEdgeCut = face.Edge(numEd)->Id();
				GenericEdge* edge = edges[idEdgeCut];
				if(edge->HasChilds())
				{
					pointsIdCutEdge[0] = edge->Point(0)->Id();
					pointsIdCutEdge[1] = edge->Point(1)->Id();
					numberEdgesChild = edge->NumberOfChilds();
					break;
				}
			}
		}
		else
		{
			GenericEdge* edge = edges[idEdge];
			pointsIdCutEdge[0] = edge->Point(0)->Id();
			pointsIdCutEdge[1] = edge->Point(1)->Id();
			numberEdgesChild = edge->NumberOfChilds();
		}

		for(unsigned int numPoint = 0; numPoint < numberEdges; numPoint++)
		{
			unsigned int idPoint = face.Point(numPoint)->Id();
			unsigned int idPointNext = face.Point((numPoint+1)%numberEdges)->Id();
			if(idPoint == pointsIdCutEdge[1] && idPointNext == pointsIdCutEdge[0])
			{
				pointsIdCutEdge[0] = idPoint;
				pointsIdCutEdge[1] = idPointNext;
				break;
			}
		}

		unsigned int numEdgesFaceChild = numberEdges - 1 + numberEdgesChild;
		faceChild.InitializeEdges(numEdgesFaceChild);
		faceChild.InitializePoints(numEdgesFaceChild);
		faceChild.AllocateFaces(numEdgesFaceChild);
		faceChild.AllocateCells(2);
		AddFace(&faceChild);

		newIdPoint.reserve(numEdgesFaceChild-1);

		//CICLO SUI LATI DELLA FACCIA
		for(unsigned int numEd = 0; numEd < numberEdges; numEd++)
		{
			unsigned int idEdge = face.Edge(numEd)->Id();
			GenericEdge* edge = edges[idEdge];
			//SE NON HA FIGLI IL LATO:
			//1) AGGIUNGO IL LATO ALLA FACCIA
			//2) AGGIORNO LA FACCIA DEL LATO PADRE CON LA FACCIA FIGLIO
			if(!edge->HasChilds())
			{
				faceChild.AddEdge(edge);
				for(unsigned int neigFace = 0; neigFace < edge->NumberOfFaces(); neigFace++)
					if(edge->Face(neigFace) != NULL)
						if(edge->Face(neigFace)->Id() == face.Id())
							edge->InsertFace(&faceChild, neigFace);
			}
			//SE HA FIGLI IL LATO:
			//1) SALVO GLI ID DEI PUNTI DEL LATO PADRE CHE E' STATO TAGLIATO
			//2) AGGIUNGO ALLA FACCIA I FIGLI DEL LATO
			//3) AI FIGLI DEL LATO AGGIORNO LA FACCIA PADRE CON LA FACCIA FIGLIO
			//4) TROVO L'ID DEL NUOVO PUNTO DA AGGIUNGERE
			//5) AGGIORNO GLI ID DEI PUNTI DEL LATO PADRE CON QUELLI DEL FIGLIO TEMPORANEAMENTE
			//   PER TROVARE L'EVENTUALE NUOVO ID DEL PUNTO SE NE SONO DUE
			//6) SALVO GLI ID DEI PUNTI DEL LATO PADRE CHE E' STATO TAGLIATO
			else
			{
				bool inversePushBack = false;
				unsigned int idEdge = edge->Child(0)->Id();
				GenericEdge* edgeChild = edges[idEdge];
				if(edgeChild->Point(0)->Id() == pointsIdCutEdge[1])
					inversePushBack = true;

				unsigned int counter = 0;
				if(inversePushBack)
				{
					for(int numChild = edge->NumberOfChilds() - 1; numChild >= 0 ; numChild--)
					{
						unsigned int idEdge = edge->Child(numChild)->Id();
						GenericEdge* edgeChild = edges[idEdge];
						faceChild.AddEdge(edgeChild);

						for(unsigned int neigFace = 0; neigFace < edge->NumberOfFaces(); neigFace++)
							if(edge->Face(neigFace) != NULL)
								if(edge->Face(neigFace)->Id() == face.Id())
									edgeChild->InsertFace(&faceChild, neigFace);

						if(counter < edge->NumberOfChilds() - 1)
						{
							if(edgeChild->Point(0)->Id() != pointsIdCutEdge[0] && edgeChild->Point(0)->Id() != pointsIdCutEdge[1])
								newIdPoint.push_back(edgeChild->Point(0)->Id());
							else
								newIdPoint.push_back(edgeChild->Point(1)->Id());
						}
						pointsIdCutEdge[0] = edgeChild->Point(0)->Id();
						pointsIdCutEdge[1] = edgeChild->Point(1)->Id();
						counter++;
					}
					pointsIdCutEdge[0] = edge->Point(0)->Id();
					pointsIdCutEdge[1] = edge->Point(1)->Id();
				}
				else
				{
					for(unsigned int numChild =  0; numChild < edge->NumberOfChilds() ; numChild++)
					{
						unsigned int idEdge = edge->Child(numChild)->Id();
						GenericEdge* edgeChild = edges[idEdge];
						faceChild.AddEdge(edgeChild);

						for(unsigned int neigFace = 0; neigFace < edge->NumberOfFaces(); neigFace++)
							if(edge->Face(neigFace) != NULL)
								if(edge->Face(neigFace)->Id() == face.Id())
									edgeChild->InsertFace(&faceChild, neigFace);

						if(counter < edge->NumberOfChilds() - 1)
						{
							if(edgeChild->Point(0)->Id() != pointsIdCutEdge[0] && edgeChild->Point(0)->Id() != pointsIdCutEdge[1])
								newIdPoint.push_back(edgeChild->Point(0)->Id());
							else
								newIdPoint.push_back(edgeChild->Point(1)->Id());
						}
						pointsIdCutEdge[0] = edgeChild->Point(0)->Id();
						pointsIdCutEdge[1] = edgeChild->Point(1)->Id();
						counter++;
					}
				}
				pointsIdCutEdge[0] = edge->Point(1)->Id();
				pointsIdCutEdge[1] = edge->Point(0)->Id();
			}
		}


		for(unsigned int numPoint = 0; numPoint < numberEdges; numPoint++)
		{
			unsigned int idPoint = face.Point(numPoint)->Id();
			unsigned int idPointNext = face.Point((numPoint+1)%numberEdges)->Id();
			GenericPoint* point = points[idPoint];
			faceChild.AddPoint(point);
			for(unsigned int facePosition = 0; facePosition < point->NumberOfFaces(); facePosition++)
			{
				//Tolgo il padre dai punti ed inserisco il figlio
				if(face.Id() == point->Face(facePosition)->Id())
					point->InsertFace(&faceChild, facePosition);
			}
			if((idPoint == pointsIdCutEdge[0] && idPointNext == pointsIdCutEdge[1]) || (idPoint == pointsIdCutEdge[1] && idPointNext == pointsIdCutEdge[0]))
			{
				for(unsigned int newId = 0; newId < newIdPoint.size() ; newId++)
				{
					GenericPoint* newPoint = points[newIdPoint[newId]];
					faceChild.AddPoint(newPoint);
					for(unsigned int facePosition = 0; facePosition < newPoint->NumberOfFaces(); facePosition++)
					{
						//Tolgo il padre dai punti ed inserisco il figlio
						if(face.Id() == newPoint->Face(facePosition)->Id())
							newPoint->InsertFace(&faceChild,facePosition);
					}
				}
			}
		}

		for(unsigned int numEd = 0; numEd < numEdgesFaceChild; numEd++)
		{
			const GenericEdge& edge = *faceChild.Edge(numEd);
			for(unsigned int numNeig = 0; numNeig < edge.NumberOfFaces(); numNeig++)
			{
				if(edge.Face(numNeig) != NULL)
				{
					GenericFace& faceNeigh = *faces[edge.Face(numNeig)->Id()];
					if(faceNeigh.Id() != faceChild.Id())
						faceChild.InsertFace(&faceNeigh, numEd);
					for(unsigned int neig = 0; neig < faceNeigh.NumberOfFaces(); neig++)
						if(faceNeigh.Face(neig) != NULL && face.Id() == faceNeigh.Face(neig)->Id())
							faceNeigh.InsertFace(&faceChild, neig);
				}
			}
		}

		for(unsigned int numCell = 0; numCell < 2; numCell++)
		{
			if(face.Cell(numCell) != NULL)
			{
				faceChild.InsertCell(face.Cell(numCell), numCell);
			}
		}
		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericMesh::CreateFaceChild(GenericFace& face, GenericFace& faceFather, const list<unsigned int>& idEdgesFace, const list<unsigned int>& idPointFace, const bool& property)
	{
		//CONTROLLARE CHE AI LATI CI SIA TUTTO
		faceFather.AddChild(&face);
		faceFather.SetState(false);

		face.SetFather(&faceFather);
		unsigned int numEdgesFace = idEdgesFace.size();
		face.InitializePoints(numEdgesFace);
		face.AllocateCells(2);
		face.AllocateEdges(numEdgesFace);
		face.AllocateFaces(numEdgesFace);
		if(property)
			face.InheritPropertiesByFather();

		//CICLO SUGLI ID DEI PUNTI DA AGGIUNGERE ALLA PAGINA
		//1) Aggiungo il punto alla faccia
		//2) Aggiorno la faccia padre con la faccia figlio ed inserisco update = true
		//3) Se Update rimane false significa che il punto è nuovo e non ha la faccia padre e bisogna
		//   aggiungere direttamente il figlio
		for(list<unsigned int>::const_iterator iteratorId = idPointFace.begin(); iteratorId != idPointFace.end(); iteratorId++)
		{
			GenericPoint* point = points[*iteratorId];
			face.AddPoint(point);
			bool update = false;
			for(unsigned int facePosition = 0; facePosition < point->NumberOfFaces(); facePosition++)
			{
				//Tolgo il padre dai punti ed inserisco il figlio
				if(faceFather.Id() == point->Face(facePosition)->Id())
				{
					update = true;
					point->InsertFace(&face, facePosition);
					break;
				}
			}
			if(!update)
				point->AddFace(&face);
		}

		unsigned int numEdge = 0;
		for(list<unsigned int>::const_iterator iteratorId = idEdgesFace.begin(); iteratorId != idEdgesFace.end(); iteratorId++)
		{
			GenericEdge* edge = edges[*iteratorId];
			face.InsertEdge(edge, numEdge);

			//CICLO SULLE FACCE VICINE AL LATO:
			//1) AGGIORNARE LA FACCIA DEL LATO CON LA FACCUA FIGLIO
			//2) AGGIUNGO LA FACCIA VICINA ALLA FACCIA CHE STO CREANDO
			for(unsigned int neigFaceEdge = 0;  neigFaceEdge < edge->NumberOfFaces(); neigFaceEdge++)
			{
				if(edge->Face(neigFaceEdge) != NULL)
					if(edge->Face(neigFaceEdge)->Id() == faceFather.Id())
						edge->InsertFace(&face, neigFaceEdge);
			}

			for(unsigned int neigFaceEdge = 0;  neigFaceEdge < edge->NumberOfFaces(); neigFaceEdge++)
			{
				if((edge->Face(neigFaceEdge) != NULL))
				{
					if((edge->Face(neigFaceEdge)->Id() != faceFather.Id()) && (edge->Face(neigFaceEdge)->Id() != face.Id()))
					{
						face.InsertFace(faces[edge->Face(neigFaceEdge)->Id()], numEdge++);
						GenericFace& faceNeigh = *faces[edge->Face(neigFaceEdge)->Id()];
						for(unsigned int numNeigh = 0; numNeigh < faceNeigh.NumberOfFaces(); numNeigh++)
							if(faceNeigh.Face(numNeigh) != NULL && faceNeigh.Face(numNeigh)->Id() == faceFather.Id())
								faceNeigh.InsertFace(&face, numNeigh);
					}
				}
				else
					numEdge++;
			}
		}

		if(faceFather.Cell(0) != NULL)
			face.InsertCell(faceFather.Cell(0), 0);

		if(faceFather.Cell(1) != NULL)
			face.InsertCell(faceFather.Cell(1), 1);
		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericMesh::ActivateFatherNodes()
	{
		for(unsigned int numEdge = 0; numEdge < edges.size(); numEdge++)
		{
			GenericEdge& edge = *edges[numEdge];
			if(edge.HasFather())
			{
				edge.SetState(false);
				const unsigned int idFather = edge.Father()->Id();
				GenericEdge& edgeFather = *edges[idFather];
				edgeFather.SetState(true);
			}
		}
		for(unsigned int numFaces = 0; numFaces < faces.size(); numFaces++)
		{
			GenericFace& face = *faces[numFaces];
			if(face.HasFather())
			{
				face.SetState(false);
				const unsigned int idFather = face.Father()->Id();
				GenericFace& faceFather = *faces[idFather];
				faceFather.SetState(true);
			}
		}
		for(unsigned int numCell = 0; numCell < cells.size(); numCell++)
		{
			GenericCell& cell = *cells[numCell];
			if(cell.HasFather())
			{
				cell.SetState(false);
				const unsigned int idFather = cell.Father()->Id();
				GenericCell& cellFather = *cells[idFather];
				cellFather.SetState(true);
			}
		}
		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericMesh::ActivateChildrenNodes()
	{
		for(unsigned int numEdge = 0; numEdge < edges.size(); numEdge++)
		{
			GenericEdge& edge = *edges[numEdge];
			if(edge.HasChilds())
			{
				edge.SetState(false);
				for(unsigned int numChild = 0; numChild < edge.NumberOfChilds(); numChild++)
				{
					const unsigned int idFather = edge.Child(numChild)->Id();
					GenericEdge& edgeChild = *edges[idFather];
					edgeChild.SetState(true);
				}
			}
		}
		for(unsigned int numFaces = 0; numFaces < faces.size(); numFaces++)
		{
			GenericFace& face = *faces[numFaces];
			if(face.HasChilds())
			{
				face.SetState(false);
				for(unsigned int numChild = 0; numChild < face.NumberOfChilds(); numChild++)
				{
					const unsigned int idFather = face.Child(numChild)->Id();
					GenericFace& faceChild = *faces[idFather];
					faceChild.SetState(true);
				}
			}
		}
		for(unsigned int numCell = 0; numCell < cells.size(); numCell++)
		{
			GenericCell& cell = *cells[numCell];
			if(cell.HasChilds())
			{
				cell.SetState(false);
				for(unsigned int numChild = 0; numChild < cell.NumberOfChilds(); numChild++)
				{
					const unsigned int idFather = cell.Child(numChild)->Id();
					GenericCell& cellChild = *cells[idFather];
					cellChild.SetState(true);
				}
			}
		}
		return Output::Success;
	}
	// ***************************************************************************
	const Output::ExitCodes GenericMesh::CleanInactiveTreeNode()
	{
		vector<GenericPoint*> pointsTemp;
		vector<GenericEdge*> edgesTemp;
		vector<GenericFace*> facesTemp;
		vector<GenericCell*> cellsTemp;

		pointsTemp.reserve(NumberOfPoints());
		edgesTemp.reserve(NumberOfEdges());
		facesTemp.reserve(NumberOfFaces());
		cellsTemp.reserve(NumberOfCells());

		unsigned int pointId = 0;
		unsigned int edgeId = 0;
		unsigned int faceId = 0;
		unsigned int cellId = 0;

		for (vector<GenericPoint*>::iterator pointPtr = points.begin(); pointPtr != points.end(); pointPtr++)
		{
			GenericPoint& point = **pointPtr;
			if(point.IsActive())
			{
				point.SetId(pointId++);
				pointsTemp.push_back(*pointPtr);
				for(int numEdge = (int)point.NumberOfEdges() - 1; numEdge >= 0; numEdge--)
					if(!point.Edge(numEdge)->IsActive())
						point.EraseEdge(numEdge);
				for(int numFace = (int)point.NumberOfFaces() - 1; numFace >= 0; numFace--)
					if(!point.Face(numFace)->IsActive())
						point.EraseFace(numFace);
				for(int numCell = (int)point.NumberOfCells() - 1; numCell >= 0 ; numCell--)
					if(!point.Cell(numCell)->IsActive())
						point.EraseCell(numCell);
			}
			else
				delete &point;
		}

		for (vector<GenericEdge*>::iterator edgePtr = edges.begin(); edgePtr != edges.end(); edgePtr++)
		{
			GenericEdge& edge = **edgePtr;
			if(edge.IsActive())
			{
				edge.SetId(edgeId++);
				edgesTemp.push_back(&edge);

				for(int numPnt = (int)edge.NumberOfPoints() - 1; numPnt >= 0; numPnt--)
					if(!edge.Point(numPnt)->IsActive())
						edge.ErasePoint(numPnt);
				for(int numEdge = (int)edge.NumberOfEdges() - 1; numEdge >= 0; numEdge--)
					if(!edge.Edge(numEdge)->IsActive())
						edge.EraseEdge(numEdge);
				for(int numFace = (int)edge.NumberOfFaces() - 1; numFace >= 0; numFace--)
					if(!edge.Face(numFace)->IsActive())
						edge.EraseFace(numFace);
				for(int numCell = (int)edge.NumberOfCells() - 1; numCell >= 0 ; numCell--)
					if((edge.Cell(numCell) != NULL) && !edge.Cell(numCell)->IsActive())
						edge.EraseCell(numCell);

				edge.SetFather(NULL);
			}
			else
				delete &edge;
		}

		for (vector<GenericFace*>::iterator facePtr = faces.begin(); facePtr != faces.end(); facePtr++)
		{
			GenericFace& face = **facePtr;
			if(face.IsActive())
			{
				face.SetId(faceId++);
				facesTemp.push_back(*facePtr);

				for(int numPnt = (int)face.NumberOfPoints() - 1; numPnt >= 0; numPnt--)
					if(!face.Point(numPnt)->IsActive())
						face.ErasePoint(numPnt);
				for(int numEdge = (int)face.NumberOfEdges() - 1; numEdge >= 0; numEdge--)
					if(!face.Edge(numEdge)->IsActive())
						face.EraseEdge(numEdge);
				for(int numFace = (int)face.NumberOfFaces() - 1; numFace >= 0; numFace--)
					if(!face.Face(numFace)->IsActive())
						face.EraseFace(numFace);
				for(int numCell = (int)face.NumberOfCells() - 1; numCell >= 0 ; numCell--)
					if((face.Cell(numCell) != NULL) && !face.Cell(numCell)->IsActive())
						face.EraseCell(numCell);

				face.SetFather(NULL);
			}
			else
				delete &face;
		}

		for (vector<GenericCell*>::iterator cellPtr = cells.begin(); cellPtr != cells.end(); cellPtr++)
		{
			GenericCell& cell = **cellPtr;
			if(cell.IsActive())
			{
				cell.SetId(cellId++);
				cellsTemp.push_back(*cellPtr);

				for(int numPnt = (int)cell.NumberOfPoints() - 1; numPnt >= 0; numPnt--)
					if(!cell.Point(numPnt)->IsActive())
						cell.ErasePoint(numPnt);
				for(int numEdge = (int)cell.NumberOfEdges() - 1; numEdge >= 0; numEdge--)
					if(!cell.Edge(numEdge)->IsActive())
						cell.EraseEdge(numEdge);
				for(int numFace = (int)cell.NumberOfFaces() - 1; numFace >= 0; numFace--)
					if(!cell.Face(numFace)->IsActive())
						cell.EraseFace(numFace);
				for(int numCell = (int)cell.NumberOfCells() - 1; numCell >= 0 ; numCell--)
					if((cell.Cell(numCell) != NULL) && !cell.Cell(numCell)->IsActive())
						cell.EraseCell(numCell);

				cell.SetFather(NULL);
			}
			else
				delete &cell;
		}

		unsigned int numPoints = pointsTemp.size();
		unsigned int numEdges = edgesTemp.size();
		unsigned int numFaces = facesTemp.size();
		unsigned int numCells = cellsTemp.size();

		points.resize(numPoints);
		edges.resize(numEdges);
		faces.resize(numFaces);
		cells.resize(numCells);

		for(unsigned int pnt = 0 ; pnt < numPoints; pnt++)
			points[pnt] = pointsTemp[pnt];

		for(unsigned int edg = 0 ; edg < numEdges; edg++)
			edges[edg] = edgesTemp[edg];

		for(unsigned int fac = 0 ; fac < numFaces; fac++)
			faces[fac] = facesTemp[fac];

		for(unsigned int cel = 0 ; cel < numCells; cel++)
			cells[cel] = cellsTemp[cel];

		return Output::Success;
	}
	// ***************************************************************************
	GenericMeshImportInterface::GenericMeshImportInterface()
	{
		maximumCellSize = 0.0;
		minimumNumberOfCells = 0;
	}
	GenericMeshImportInterface::~GenericMeshImportInterface()
	{
	}
	// ***************************************************************************
	void GenericMeshImportInterface::SetBoundaryConditions(const vector<int>& _vertexMarkers, const vector<int>& _edgeMarkers, const vector<int>& _faceMarkers)
	{
		if (_vertexMarkers.size() > 0)
		{
			vertexMarkers.resize(_vertexMarkers.size());
			memcpy(&vertexMarkers[0], &_vertexMarkers[0], _vertexMarkers.size() * sizeof(int));
		}

		if (_edgeMarkers.size() > 0)
		{
			edgeMarkers.resize(_edgeMarkers.size());
			memcpy(&edgeMarkers[0], &_edgeMarkers[0], _edgeMarkers.size() * sizeof(int));
		}

		if (_faceMarkers.size() > 0)
		{
			faceMarkers.resize(_faceMarkers.size());
			memcpy(&faceMarkers[0], &_faceMarkers[0], _faceMarkers.size() * sizeof(int));
		}
	}
	// ***************************************************************************
}
