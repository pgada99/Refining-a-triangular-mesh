#include <math.h>
#include <iomanip>
#include <sys/stat.h>

#include "MeshImport_Triangle.hpp"
#include "Output.hpp"

using namespace MainApplication;

namespace GeDiM
{
	// ***************************************************************************
	MeshImport_Triangle::MeshImport_Triangle() : GenericMeshImportInterface()
	{
		inputMeshPointer = NULL;
		outputMeshPointer = NULL;
		triangleOptions = "-Qzpqnea";
	}
	MeshImport_Triangle::~MeshImport_Triangle()
	{
		if (inputMeshPointer != NULL)
		{
			delete[] inputMeshPointer->pointlist; inputMeshPointer->pointlist = NULL;
			delete[] inputMeshPointer->pointmarkerlist; inputMeshPointer->pointmarkerlist = NULL;
			delete[] inputMeshPointer->segmentlist; inputMeshPointer->segmentlist = NULL;
			delete[] inputMeshPointer->segmentmarkerlist; inputMeshPointer->segmentmarkerlist = NULL;
		}

		if (outputMeshPointer != NULL)
		{
			free(outputMeshPointer->pointlist);
			trifree( outputMeshPointer->pointmarkerlist );
			trifree( outputMeshPointer->segmentlist );
			trifree( outputMeshPointer->segmentmarkerlist );
			free( outputMeshPointer->edgelist);
			free( outputMeshPointer->edgemarkerlist);
			free( outputMeshPointer->trianglearealist);
			free( outputMeshPointer->neighborlist);
		}

		delete inputMeshPointer; inputMeshPointer = NULL;
		delete outputMeshPointer; outputMeshPointer = NULL;
	}
	// ***************************************************************************
	Output::ExitCodes MeshImport_Triangle::CreateTriangleInput(const GenericDomain& domain)
	{
		delete inputMeshPointer; inputMeshPointer = NULL;
		inputMeshPointer = new triangulateio();

		const GenericDomain2D& domain2D = dynamic_cast<const GenericDomain2D&>(domain); //cast da GenericDomain (classe base) a GenericDomain2D (classe derivata)

		const unsigned int& numberOfVertices = domain2D.TotalNumberVertices();
		const unsigned int& numberOfEdges = domain2D.TotalNumberEdges();

		if (numberOfVertices == 0 || numberOfEdges == 0)
		{
			Output::PrintErrorMessage("Wrong initialization of the domain %d, no vertices or edges", false, domain2D.GlobalId());
			return Output::GenericError;
		}

		inputMeshPointer->pointlist = new double[2 * numberOfVertices];
		inputMeshPointer->pointattributelist = NULL;
		inputMeshPointer->pointmarkerlist = new int[numberOfVertices];
		inputMeshPointer->numberofpoints = numberOfVertices;
		inputMeshPointer->numberofpointattributes = 0;
		inputMeshPointer->numberofsegments = numberOfVertices;
		inputMeshPointer->trianglelist = NULL;
		inputMeshPointer->triangleattributelist = NULL;
		inputMeshPointer->trianglearealist = NULL;
		inputMeshPointer->neighborlist = NULL;
		inputMeshPointer->numberoftriangles = 0;
		inputMeshPointer->numberofcorners = 0;
		inputMeshPointer->numberoftriangleattributes = 0;
		inputMeshPointer->segmentlist = new int[2 * numberOfVertices];
		inputMeshPointer->segmentmarkerlist = new int[numberOfVertices];
		inputMeshPointer->holelist = NULL;
		inputMeshPointer->numberofholes = 0;
		inputMeshPointer->regionlist = NULL;
		inputMeshPointer->numberofregions = 0;
		inputMeshPointer->edgelist = NULL;
		inputMeshPointer->edgemarkerlist = NULL;
		inputMeshPointer->normlist = NULL;
		inputMeshPointer->numberofedges = 0;

		double* point_list = inputMeshPointer->pointlist;
		int* point_markerlist = inputMeshPointer->pointmarkerlist;
		int* segment_list = inputMeshPointer->segmentlist;
		int* segment_markerlist = inputMeshPointer->segmentmarkerlist;

		for (unsigned int j = 0; j < numberOfVertices; j++)
		{
			point_list[2 * j] = domain2D.RotatedVertex(j)(0);
			point_list[2 * j + 1] = domain2D.RotatedVertex(j)(1);

			point_markerlist[j] = 2;
		}

		for (unsigned int j = 0; j < numberOfEdges; j++)
		{
			segment_list[2 * j] = domain2D.EdgeOriginIndex(j);
			segment_list[2 * j + 1] = domain2D.EdgeEndIndex(j);

			segment_markerlist[j]= 2;
		}

    if (!vertexMarkers.empty())
			memcpy(point_markerlist, vertexMarkers.data(), numberOfVertices * sizeof(int));
		if (!edgeMarkers.empty())
			memcpy(segment_markerlist, edgeMarkers.data(), numberOfEdges * sizeof(int));

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes MeshImport_Triangle::CreateTriangleOutput(const GenericDomain& domain)
	{
		if (minimumNumberOfCells == 0 && maximumCellSize <= 0)
		{
			Output::PrintErrorMessage("Wrong initialization of the minimumNumberOfCells or minimumCellSize", false);
			return Output::GenericError;
		}

		if (inputMeshPointer == NULL)
		{
			Output::PrintErrorMessage("No Triangle input in domain %d", false, domain.GlobalId());
			return Output::GenericError;
		}

		const GenericDomain2D& domain2D = dynamic_cast<const GenericDomain2D&>(domain);

		if (minimumNumberOfCells > 0 && domain2D.Area() == 0)
		{
			Output::PrintErrorMessage("Wrong initialization of the domain %d", false, domain2D.GlobalId());
			return Output::GenericError;
		}

		const double& domainArea = domain2D.Area();
		double cellArea = minimumNumberOfCells == 0 ? maximumCellSize : domainArea / (double)minimumNumberOfCells;

		ostringstream options;
		options.precision(16);
		options<< triangleOptions;
		options<< cellArea;
		size_t sizeOptions = options.str().size();
		char* optionPointer = new char[sizeOptions + 1];
		options.str().copy(optionPointer, sizeOptions);
		optionPointer[sizeOptions] = '\0';

		delete outputMeshPointer; outputMeshPointer = NULL;
		outputMeshPointer = new triangulateio();

		triangulate(optionPointer, inputMeshPointer, outputMeshPointer, (struct triangulateio*)NULL);

		delete[] optionPointer;

		return Output::Success;
	}
	// ***************************************************************************
	Output::ExitCodes MeshImport_Triangle::CreateMesh(const GenericDomain& domain, GenericMesh& mesh) const
	{
		/// <ul>

		if (outputMeshPointer == NULL)
		{
			Output::PrintErrorMessage("No Triangle ouput in domain %d", false, domain.GlobalId());
			return Output::GenericError;
		}

		const GenericDomain2D& domain2D = dynamic_cast<const GenericDomain2D&>(domain);

		/// <li>	Fill mesh structures
		unsigned int numberOfCellMesh = outputMeshPointer->numberoftriangles;
		unsigned int numberOfEgdesMesh = outputMeshPointer->numberofedges;
		unsigned int numberOfPointsMesh = outputMeshPointer->numberofpoints;

		mesh.InitializeCells(numberOfCellMesh);
		mesh.InitializePoints(numberOfPointsMesh);
		mesh.InitializeEdges(numberOfEgdesMesh);
		vector<GenericCell*> cells(numberOfCellMesh);
		vector<GenericEdge*> edges(numberOfEgdesMesh);
		vector<GenericPoint*> points(numberOfPointsMesh);

		/// <li> Set Points
		for (unsigned int p = 0; p < numberOfPointsMesh; p++)
		{
			points[p] = mesh.CreatePoint();

			GenericPoint* point = points[p];

			Vector3d point2d(outputMeshPointer->pointlist[2 * p], outputMeshPointer->pointlist[2 * p + 1], 0.0);
			Vector3d rotatedPoint = domain2D.RotatePoint(point2d, true, true);

			point->SetCoordinates(rotatedPoint);
			point->SetMarker(outputMeshPointer->pointmarkerlist[p]);

			mesh.AddPoint(point);
		}

		/// <li> Set Edges
		for(unsigned int ed = 0; ed < numberOfEgdesMesh; ed++)
		{
			edges[ed] = mesh.CreateEdge();

      GenericEdge* edge = edges[ed];
      edge->AllocateCells(2);

			edge->SetMarker(outputMeshPointer->edgemarkerlist[ed]);

			for(int i = 0; i < 2; i++)
			{
				GenericPoint* point = points[outputMeshPointer->edgelist[2 * ed + i]];
				edge->AddPoint(point);
				point->AddEdge(edge);
			}

			mesh.AddEdge(edge);
		}

		/// <li> Set Cells
    for (unsigned int c = 0; c < numberOfCellMesh; c++)
		{
			cells[c] = mesh.CreateCell();

			GenericCell* cell = cells[c];
			cell->AllocatePoints(3); //alloca un vettore per i 3 vertici
			cell->AllocateEdges(3); //alloca un vettore per i 3 lati
			cell->AllocateCells(3); //alloca un vettore per le 3 celle vicine (triangoli vicini)

			for (int i = 0; i < 3; i++)
			{
				GenericPoint* point = points[outputMeshPointer->trianglelist[3 * c + i]];

				cell->InsertPoint(point,i); //aggiunge i punti al triangolo
				point->AddCell(cell); //aggiunge il triangolo ai punti
			}
			mesh.AddCell(cell); // aggiunge la cella alla mesh
		}

    //Inserisce nelle celle le informazioni sulle celle vicine
		for (unsigned int c = 0; c < numberOfCellMesh; c++)
		{
			GenericCell* cell = cells[c];
			for (int i = 0; i < 3; i++)
			{
				if (outputMeshPointer->neighborlist[3 * c + i] > -1)
				{
					GenericCell* cellNeigh = cells[outputMeshPointer->neighborlist[3 * c + i]];
					cell->InsertCell(cellNeigh,(i+1)%3);
				}
			}
		}

		vector<int> edgeVerteceIds(2, -1);
		vector<int> cellVertceIds(3, -1);

		for(int ed = 0; ed < outputMeshPointer->numberofedges; ed++)
		{
			GenericEdge* edge = edges[ed];
			edgeVerteceIds[0] = outputMeshPointer->edgelist[2*ed];
			edgeVerteceIds[1] = outputMeshPointer->edgelist[2*ed + 1];

			for(int e = 0; e < outputMeshPointer->numberoftriangles; e++)
			{
				cellVertceIds[0] = outputMeshPointer->trianglelist[3*e];
				cellVertceIds[1] = outputMeshPointer->trianglelist[3*e + 1];
				cellVertceIds[2] = outputMeshPointer->trianglelist[3*e + 2];

				int k = CheckTrianglePosition(edgeVerteceIds, cellVertceIds);

				if(k != 0)
				{
					GenericCell* cell = cells[e];
					if(k>9 && k<13)
					{
						edge->InsertCell(cell,0) ; // right
						if(k==10)
							cell->InsertEdge(edge,0);
						else if(k==11)
							cell->InsertEdge(edge,1);
						else if(k==12)
							cell->InsertEdge(edge,2);
					}
					if(k>19 && k<23)
					{
						edge->InsertCell(cell,1) ; // left
						if(k==20)
							cell->InsertEdge(edge,0);
						else if(k==21)
							cell->InsertEdge(edge,1);
						else if(k==22)
							cell->InsertEdge(edge,2);
					}
				}
			}
		}

		/// <li> Check cell edges
		for (int e = 0; e < outputMeshPointer->numberoftriangles; e++)
		{
      for (int i = 0; i < 3; i++)
			{
				if (cells[e]->Edge(i) == NULL)
				{
					Output::PrintErrorMessage("Error on meshing on domain %d: cell %d miss edge number %d", false, domain.GlobalId(), e, (i + 1));
					exit(-1);
				}
			}
		}

		/// <li> Check edge points
		for (int ed = 0; ed < outputMeshPointer->numberofedges; ed++)
		{
      for (int i=0; i < 2; i++)
			{
				if (edges[ed]->Point(i) == NULL)
				{
					Output::PrintErrorMessage("Error on meshing on domain %d: edge %d miss point number %d", false, domain.GlobalId(), ed, (i + 1));
					exit(-1);
				}
			}
		}

		return Output::Success;

		/// </ul>
	}
	// Usando l'ordinamento in senso antiorario dei vertici del triangolo, trova la posizione del triangolo vicino
	int MeshImport_Triangle::CheckTrianglePosition(const vector<int>& edgePointIds,const vector<int>& cellPointIds)
	{
		if(edgePointIds[0]==cellPointIds[2] && edgePointIds[1]==cellPointIds[1])
			return 11;  //right, edge 1
		else if(edgePointIds[0]==cellPointIds[1] && edgePointIds[1]==cellPointIds[0])
			return 10;  //right, edge 0
		else if(edgePointIds[0]==cellPointIds[0] && edgePointIds[1]==cellPointIds[2])
			return 12;  //right, edge 2
		else if(edgePointIds[0]==cellPointIds[0] && edgePointIds[1]==cellPointIds[1])
			return 20; //left, edge 0
		else if(edgePointIds[0]==cellPointIds[1] && edgePointIds[1]==cellPointIds[2])
			return 21; //left, edge 1
		else if(edgePointIds[0]==cellPointIds[2] && edgePointIds[1]==cellPointIds[0])
			return 22; //left, edge 2

		return 0; // edge not adjacent to cell
	}
	// ***************************************************************************
	Output::ExitCodes MeshImport_Triangle::ExportTriangleMesh(const string& nameFolder, const string& nameFile) const
	{
		if (inputMeshPointer == NULL || outputMeshPointer == NULL)
			return Output::GenericError;

		ostringstream nameFolderStream, nameFileStream;

		nameFolderStream<< nameFolder<< "/";
		nameFolderStream<< "Triangle/";

		Output::CreateFolder(nameFolderStream.str());

		nameFileStream<< nameFolderStream.str()<< nameFile;

		const struct triangulateio& input = *inputMeshPointer;
		const struct triangulateio& triangulation = *outputMeshPointer;


		ofstream file;
		char basefilename[50];
		char filename[50];
		strcpy(basefilename, nameFileStream.str().c_str());

		file<< setprecision(16);

		strcpy (filename, basefilename);
		strcat (filename, ".poly");
		file.open(filename);
		file<< input.numberofpoints<< " "<< "2"<< " "<< "0"<< " "<< "1"<< endl;
		for(int i = 0; i < input.numberofpoints; i++)
			file<< i + 1<< " "
					<< input.pointlist[2 * i]<< " "
					<< input.pointlist[2 * i + 1]<< " "
					<< input.pointmarkerlist[i]<< endl;

		file<< input.numberofsegments<< " "<< " 0 "<< endl;

		for(int i = 0; i < input.numberofsegments; i++)
			file<< i + 1<< " "
					<< input.segmentlist[2 * i] + 1<< " "
					<< input.segmentlist[2 * i + 1] + 1<< " "
					<< input.segmentmarkerlist[i]<< endl;

		file<< "0";
		file.close();

		strcpy (filename, basefilename);
		strcat (filename, ".1.node");
		file.open(filename);
		file<< triangulation.numberofpoints<< " "<< "2"<< " "<< "0"<< " "<< "1"<< endl;
		for(int i = 0; i < triangulation.numberofpoints; i++)
			file<< i + 1<< " "<< triangulation.pointlist[2 * i]<< " "
					<< triangulation.pointlist[2 * i + 1]<< " "
					<< triangulation.pointmarkerlist[i]<< endl;
		file.close();

		strcpy (filename, basefilename);
		strcat (filename, ".1.neigh");
		file.open(filename);
		file<< triangulation.numberoftriangles<< " "<< "3"<< endl;
		for(int i = 0; i < triangulation.numberoftriangles; i++)
			file<< i + 1<< " "
					<< ((triangulation.neighborlist[3 * i] != -1) ? (triangulation.neighborlist[3 * i] + 1) : -1)<< " "
					<< ((triangulation.neighborlist[3 * i + 1] != -1) ? (triangulation.neighborlist[3 * i + 1] + 1) : -1)<< " "
					<< ((triangulation.neighborlist[3 * i + 2] != -1) ? (triangulation.neighborlist[3 * i + 2] + 1) : -1)<< endl;
		file.close();

		strcpy (filename,basefilename);
		strcat (filename,".1.edge");
		file.open(filename);
		file<< triangulation.numberofedges<< " "<< "1"<< endl;
		for(int i = 0; i < triangulation.numberofedges; i++)
			file<< i + 1<< " "<< triangulation.edgelist[2 * i] + 1
					<< " "<< triangulation.edgelist[2 * i + 1] + 1
					<< " "<< triangulation.edgemarkerlist[i]<< endl;
		file.close();

		strcpy (filename,basefilename);
		strcat (filename,".1.ele");
		file.open(filename);
		file<< triangulation.numberoftriangles<< " "<< "3"<< " "<< "0"<< endl;
		for(int i = 0; i < triangulation.numberoftriangles; i++)
			file<< i + 1<< " "
					<< triangulation.trianglelist[3 * i] + 1<< " "
					<< triangulation.trianglelist[3 * i + 1] + 1<< " "
					<< triangulation.trianglelist[3 * i + 2] + 1<< endl;
		file.close();

		strcpy (filename,basefilename);
		strcat (filename,".1.poly");
		file.open(filename);
		file<< " 0 "<< "2 "<< "0 "<< "1"<< endl;
		file<< triangulation.numberofsegments<< " 1"<< endl;
		for(int i = 0; i < triangulation.numberofsegments; i++)
		{
			file<< i + 1<< " "
					<< triangulation.segmentlist[2 * i] + 1<< " "
					<< triangulation.segmentlist[2 * i + 1] + 1<< " "
					<< triangulation.segmentmarkerlist[i] << endl;
		}

		file<< "0";
		file.close();

		return Output::Success;
	}
	// ***************************************************************************
}
