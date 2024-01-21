#include "TriangleRefiner.hpp"
#include <math.h>
#include <string>
//#define DEBUG_EDGE 1
//#define DEBUG_REFINER 1
//#define BREAKPOINTS 1
//#define DEBUG_NEIGHBORHOOD 1

using namespace std;

/*
 *  NOTE
 * 	-	L'orientazione dei lati viene calcolata in base all'ordine dei punti nella cella (GenericCell) e non in base all'ordine dei punti
 * 		memorizzati nel GenericEdge poichè ciascun lato è condiviso da due triangoli e per ciascuno di essi avrà orientazione diversa.
 *
 *	TODO
 *	-	Aggiornare la vicinanza anche per i punti e per le celle (point-edge, point-cell, cell-cell)
 *
 *	DUBBI
 *	-	Ho messo la vicinanza point-edge, point-cell, cell-cell. Ha senso mettere la vicinanza edge-edge?
 */

namespace GeDiM{

	/*
	 * Funzione che scorre il vettore idCellsToRefine e chiama, per ogni triangolo,
	 * la funzione che lo taglia. Se il triangolo è inattivo non viene tagliato una
	 * seconda volta.
	 */
	Output::ExitCodes TriangleRefiner::RefineMesh(){
		unsigned int i=0;
		while(i<idCellsToRefine.size()){

			if(idCellsToRefine.at(i) <= meshPointer->NumberOfCells()){ // Se il numero della cella da rifinire è più piccolo del numero totale di celle nella mesh
				/*
				 * Verifica che il lato sia attivo per elaborare solo le foglie dell'albero
				 */
				if( (meshPointer->Cell(idCellsToRefine.at(i)))->IsActive() ){
					RefineTriangle(idCellsToRefine.at(i));
				}
			}
			i++;
		}
		idCellsToRefine.clear();
		UpdateNeighbourhood();
		//PrintNeigh();
		return Output::Success;
	}

	/*
	 * Funzione che ruota i punti ed i lati del triangolo in modo tale che il lato E0 sia il più lungo
	 * e che i lati successivi siano ordinati secondo il verso antiorario.
	 */
	void TriangleRefiner::Rotate( const unsigned int& value ){
		GenericEdge *E0, *E1, *E2, *aux; // lati del triangolo
		GenericPoint *P0, *P1, *P2, *auxp; // punti del triangolo

		GenericCell* cell=meshPointer->Cell(value);

		E0=meshPointer->Edge(cell->Edge(0)->Id());
		E1=meshPointer->Edge(cell->Edge(1)->Id());
		E2=meshPointer->Edge(cell->Edge(2)->Id());

		P0=meshPointer->Point(cell->Point(0)->Id());
		P1=meshPointer->Point(cell->Point(1)->Id());
		P2=meshPointer->Point(cell->Point(2)->Id());

#ifdef BREAKPOINTS
		cout<<"ROT BP1"<<endl;
#endif
		/*
		 * Calcolo del quadrato della lunghezza dei 3 lati: non viene calcolata la lunghezza
		 * poichè il calcolo della radice quadrata ha un costo computazionale elevato e per il
		 * confronto è sufficiente il quadrato.
		 */
		double length[3];
		length[0]=pow(E0->Point(0)->X()-E0->Point(1)->X(),2)+pow(E0->Point(0)->Y()-E0->Point(1)->Y(),2);
		length[1]=pow(E1->Point(0)->X()-E1->Point(1)->X(),2)+pow(E1->Point(0)->Y()-E1->Point(1)->Y(),2);
		length[2]=pow(E2->Point(0)->X()-E2->Point(1)->X(),2)+pow(E2->Point(0)->Y()-E2->Point(1)->Y(),2);

		/*
		 * Ricerca del lato più lungo e scambio dei lati nella cella in modo da avere il lato E0 come lato più lungo
		 */
		if(length[1]==length[2] && length[1]==length[0]){
					#ifdef DEBUG_EDGE
						cout<<"Triangolo equilatero."<<endl;
					#endif
		} else if(length[1]>=length[2] && length[1]>=length[0]){
			cell->InsertEdge(E1,0);
			cell->InsertEdge(E2,1);
			cell->InsertEdge(E0,2);

			aux=E0;
			E0=E1;
			E1=E2;
			E2=aux;
					#ifdef DEBUG_EDGE
						cout<<"Il lato 1 è il maggiore."<<endl;
						cout<<"Lato 0: "<<length[0]<<", Lato 1: "<<length[1]<<", Lato 2: "<<length[2]<<endl;
					#endif
		} else if(length[2]>=length[0] && length[2]>=length[1]){
			cell->InsertEdge(E2,0);
			cell->InsertEdge(E0,1);
			cell->InsertEdge(E1,2);

			aux=E1;
			E1=E0;
			E0=E2;
			E2=aux;
					#ifdef DEBUG_EDGE
						cout<<"Il lato 2 è il maggiore."<<endl;
						cout<<"Lato 0: "<<length[0]<<", Lato 1: "<<length[1]<<", Lato 2: "<<length[2]<<endl;
					#endif
		} else if (length[0]>=length[1] && length[0]>=length[2]){
					#ifdef DEBUG_EDGE
						cout<<"Il lato 0 è il più lungo."<<endl;
					#endif
		}

		/*
		 * Ricerca dei punti P0, P1 e P2 del triangolo in base al nuovo ordine dei lati
		 */
#ifdef BREAKPOINTS
		cout<<"ROT BP2"<<endl;
#endif
		//P0
		if( E0->Point(0)->Id() == E2->Point(0)->Id() || E0->Point(0)->Id() == E2->Point(1)->Id() ){
			P0=meshPointer->Point(E0->Point(0)->Id());
		} else if ( E0->Point(1)->Id() == E2->Point(0)->Id() || E0->Point(1)->Id() == E2->Point(1)->Id() ){
			P0=meshPointer->Point(E0->Point(1)->Id());
		}

		//P1
		if( E1->Point(0)->Id() == E0->Point(0)->Id() || E1->Point(0)->Id() == E0->Point(1)->Id() ){
			P1=meshPointer->Point(E1->Point(0)->Id());
		} else if ( E1->Point(1)->Id() == E0->Point(0)->Id() || E1->Point(1)->Id() == E0->Point(1)->Id() ){
			P1=meshPointer->Point(E1->Point(1)->Id());
		}

		//P2
		if( E2->Point(0)->Id() == E1->Point(0)->Id() || E2->Point(0)->Id() == E1->Point(1)->Id() ){
			P2=meshPointer->Point(E2->Point(0)->Id());
		} else if ( E2->Point(1)->Id() == E1->Point(0)->Id() || E2->Point(1)->Id() == E1->Point(1)->Id() ){
			P2=meshPointer->Point(E2->Point(1)->Id());
		}


		/*
		 * Calcolo del prodotto vettoriale tra i lati E0 ed E2 del triangolo per trovare l'orientazione
		 * del lato E0: se il prodotto vettoriale dà risultato negativo significa che l'orientazione
		 * attuale del lato E0 è errata, quindi si scambiano i punti P0 e P1 del triangolo.
		 */

		if( ( ( (P1->X()-P0->X())*(P2->Y()-P0->Y()) ) - ( (P1->Y()-P0->Y())*(P2->X()-P0->X()) ) ) < 0 ){
			auxp=P1;
			P1=P0;
			P0=auxp;
			cout<<"cella "<<value<<" PV neg"<<endl;
		}

		/*
		 * Memorizza i punti nella cella secondo l'ordine corretto
		 */
		cell->InsertPoint(P0,0);
		cell->InsertPoint(P1,1);
		cell->InsertPoint(P2,2);

#ifdef BREAKPOINTS
		cout<<"ROT BP3"<<endl;
#endif
		return;
	}


	/*
	 * Funzione che taglia il triangolo
	 */
	  Output::ExitCodes TriangleRefiner::RefineTriangle(const unsigned int& value){
		#ifdef DEBUG_REFINER
		  cout<<"\n\nCella da rifinire: "<<value<<endl;
		#endif
		  Rotate(value);

		  GenericCell* cell = meshPointer->Cell(value); // Puntatore alla cella di indice value (DA DIVIDERE)
		  GenericEdge* E0 = meshPointer->Edge(cell->Edge(0)->Id()); // Puntatore al lato 0 della cella-> prende l'id del lato (È UN PUNTATORE COSTANTE, NON POSSO LAVORARCI DIRETTAMENTE)
		  GenericEdge* E1 = meshPointer->Edge(cell->Edge(1)->Id());
		  GenericEdge* E2 = meshPointer->Edge(cell->Edge(2)->Id());
		  GenericEdge* Em = meshPointer->CreateEdge(); // Nuovo lato che verrà generato
		  meshPointer->AddEdge(Em); // Aggiunge il puntatore alla Mesh: le informazioni sul lato verranno calcolate e aggiunte in seguito

		  GenericPoint* P0 = meshPointer->Point(cell->Point(0)->Id()); // Puntatori ai vertici della cella
		  GenericPoint* P1 = meshPointer->Point(cell->Point(1)->Id());
		  GenericPoint* P2 = meshPointer->Point(cell->Point(2)->Id());
#ifdef DEBUG_NEIGHBORHOOD
		  cout<<"\tCoordinate: P0 ("<<P0->X()<<", "<<P0->Y()<<") - P1 ("<<P1->X()<<", "<<P1->Y()<<") - P2 ("<<P2->X()<<", "<<P2->Y()<<")"<<endl;
		  cout<<"\tLati di "<<cell->Id()<<": "<<E0->Id()<<", "<<E1->Id()<<" e "<<E2->Id()<<endl;
#endif

		  GenericPoint* Pm;

		  GenericEdge* Ei;
		  GenericEdge* Ee;
#ifdef BREAKPOINTS
		  cout<<"RT BP1"<<endl;
#endif

#ifdef DEBUG_NEIGHBORHOOD
		  string right;
		  string left;
		  if(E0->HasRightCell()){
			  right=to_string(E0->RightCell()->Id());
		  } else {
			  right="void";
		  }
		  if(E0->HasLeftCell()){
			  left=to_string(E0->LeftCell()->Id());
		  } else {
			  left="void";
		  }
		  cout<<"\t\tVicini di E0 [id="<<E0->Id()<<"]: "<<right<<" e "<<left<<endl;

		  if(E1->HasRightCell()){
			  right=to_string(E1->RightCell()->Id());
		  } else {
			  right="void";
		  }
		  if(E1->HasLeftCell()){
			  left=to_string(E1->LeftCell()->Id());
		  } else {
			  left="void";
		  }
		  cout<<"\t\tVicini di E1 [id="<<E1->Id()<<"]: "<<right<<" e "<<left<<endl;

		  if(E2->HasRightCell()){
			  right=to_string(E2->RightCell()->Id());
		  } else {
			  right="void";
		  }
		  if(E2->HasLeftCell()){
			  left=to_string(E2->LeftCell()->Id());
		  } else {
			  left="void";
		  }
		cout<<"\t\tVicini di E2 [id="<<E2->Id()<<"]: "<<right<<" e "<<left<<endl;
#endif

		  if( E0->IsActive() ){ 		// Verifica che il lato considerato sia attivo (non già rifinito)
				#ifdef DEBUG_NEIGHBORHOOD
							  cout<<"LATO E0 ATTIVO: TAGLIO LATO E CELLA"<<endl;
				#endif
			  Ei = meshPointer->CreateEdge(); // Crea il nuovo lato da P0 al Punto medio Pm
			  meshPointer->AddEdge(Ei);
			  Ee = meshPointer->CreateEdge(); // Crea il nuovo lato da P0 al Punto medio Pm
			  meshPointer->AddEdge(Ee);

			  Pm = meshPointer->CreatePoint(); // Crea il puntatore al punto medio del lato che verrà tagliato
			  meshPointer->AddPoint(Pm); // Aggiunge il punto al vettore di punti della Mesh
			  Pm->SetCoordinates( 0.5*(P0->X()+P1->X()), 0.5*(P0->Y()+P1->Y()) ); // Calcola le coordinate del punto medio
			  Pm->SetMarker(E0->Marker()); // Imposta il Marker del punto prendendolo dal marker del lato

			  /*
			   * Crea i due nuovi Edge che sostituiranno l'Edge E0 tagliato
			   * Ei: prima metà del lato tagliato, Ee: seconda metà del lato tagliato
			   */
			  Ei->AddPoint(P0);
			  Ei->AddPoint(Pm);

			  Ee->AddPoint(Pm);
			  Ee->AddPoint(P1);

			  /*
			   * Aggiunge i punti Pm e P2 ad Em
			   */
			  Em->AddPoint(Pm);
			  Em->AddPoint(P2);

			  /*
			   *  Gestisce l'inserimento dei figli nell'albero
			   */
			  E0->SetState(false); // E0 viene disattivato poichè è stato diviso
			  E0->InitializeChilds(2);
			  E0->AddChild(Ei); // Ei è il primo figlio
			  E0->AddChild(Ee); // Ee è il secondo figlio
			  Ei->SetFather(E0);
			  Ee->SetFather(E0);
			  Ei->SetMarker(E0->Marker());
			  Ee->SetMarker(E0->Marker());

			  Ee->InitializeCells(2);
			  Ei->InitializeCells(2);

#ifdef BREAKPOINTS
	  cout<<"RT BP2"<<endl;
#endif
			  /*
			   * Dopo aver tagliato la cella, aggiunge la vicina alle celle da tagliare
			   */
			  if(E0->HasRightCell() && E0->RightCell()->Id()==cell->Id() && E0->HasLeftCell()){
					AddCellId(E0->LeftCell()->Id());
					#ifdef DEBUG_REFINER
						cout<<"\tAggiungo la cella "<<E0->LeftCell()->Id()<<" (cella vicina sx)"<<endl;
					#endif
			  }

			  if(E0->HasLeftCell() && E0->LeftCell()->Id()==cell->Id() && E0->HasRightCell()){
					AddCellId(E0->RightCell()->Id());
					#ifdef DEBUG_REFINER
						 cout<<"\tAggiungo la cella "<<E0->RightCell()->Id()<<" (cella vicina dx)"<<endl;
					#endif
			  }

				#ifdef BREAKPOINTS
							  cout<<"RT BP3"<<endl;
				#endif

		  } else { // Se il lato E0 è già stato tagliato -> taglio solo la cella
				#ifdef DEBUG_NEIGHBORHOOD
						cout<<"LATO E0 INATTIVO: TAGLIO SOLO LA CELLA"<<endl;
				#endif

			  /*
			   * Nel caso in cui tutti e tre i lati siano già inattivi, il triangolo viene
			   * tagliato in 4 triangoli simili e si termina la funzione.
			   */
			  if(!E1->IsActive() && !E2->IsActive()){
				  RefineUniformly(value);
				  #ifdef DEBUG_REFINER
				  	  cout<<"Refine uniformly (cella "<<cell->Id()<<")"<<endl;
				  #endif
				  return Output::Success;
			  }

#ifdef BREAKPOINTS
			  cout<<"RT BP4"<<endl;
#endif
			  /*
			   * Trova il punto medio del lato E0 cercando il punto comune tra i due lati
			   * figli di E0. [Point0<-Child0]
			   */
			  GenericPoint *P0C0=meshPointer->Point( ( (GenericEdge*)(E0->Child(0)) ) ->Point(0)->Id());
			  GenericPoint *P1C0=meshPointer->Point( ( (GenericEdge*)(E0->Child(0)) ) ->Point(1)->Id());
			  GenericPoint *P0C1=meshPointer->Point( ( (GenericEdge*)(E0->Child(1)) ) ->Point(0)->Id());
			  GenericPoint *P1C1=meshPointer->Point( ( (GenericEdge*)(E0->Child(1)) ) ->Point(1)->Id());

			  if( P0C0->Id() == P0C1->Id() || P0C0->Id() == P1C1->Id() ){
				  Pm=P0C0;
			  } else if ( P1C0->Id() == P0C1->Id() || P1C0->Id() == P1C1->Id() ){
				  Pm=P1C0;
			  }
			  cout<<"Pm ("<<Pm->X()<<", "<<Pm->Y()<<")"<<endl;

			  /*
			   * Aggiunge le informazioni sui punti di Em (lato che taglia a metà il triangolo)
			   */
			  Em->AddPoint(Pm);
			  Em->AddPoint(cell->Point(2));

			  /*
			   * Ei ed Ee vengono presi al contrario poichè il lato è stato tagliato guardando dal triangolo
			   * accanto e quindi i due lati sono stati inseriti in ordine inverso:
			   * Ei -> secondo figlio
			   * Ee -> primo figlio
			   */
			  Ei=meshPointer->Edge(E0->Child(1)->Id());
			  Ee=meshPointer->Edge(E0->Child(0)->Id());

#ifdef BREAKPOINTS
			  cout<<"RT BP5"<<endl;
#endif
		  } // if(E0->IsActive()) {} else {}


		  /*
		   * Crea le due nuove celle e le aggiunge alla mesh
		   */
		  GenericCell* Cd=meshPointer->CreateCell();
		  meshPointer->AddCell(Cd);
		  GenericCell* Cs=meshPointer->CreateCell();
		  meshPointer->AddCell(Cs);

		  cell->SetState(false); // Disattiva la cella che è stata tagliata
		  cell->InitializeChilds(2); // Crea e aggiunge i figli
		  cell->AddChild(Cd);
		  cell->AddChild(Cs);
		  Cd->SetFather(cell);
		  Cs->SetFather(cell);

#ifdef BREAKPOINTS
		  cout<<"RT BP6"<<endl;
#endif

		  /*
		   * Poichè come prima operazione ruotiamo la cella da rifinire ed ordiniamo i suoi lati ed i suoi punti, e poichè utilizziamo solo
		   * le informazioni della cella per effettuare le operazioni, non è necessario verificare la corrispondenza tra P0 ed E0->Point(0)
		   */

		  /*
		   * Inizializza e aggiunge i 3 punti della cella destra
		   */
		  Cd->InitializePoints(3);
		  Cd->AddPoint(P0);
		  Cd->AddPoint(Pm);
		  Cd->AddPoint(P2);

		  /*
		   * Inizializza e aggiunge i 3 punti della cella di sinistra
		   */
		  Cs->InitializePoints(3);
		  Cs->AddPoint(Pm);
		  Cs->AddPoint(P1);
		  Cs->AddPoint(P2);

		  /*
		   * Inizializza e aggiunge i 3 lati della cella di destra
		   */
		  Cd->InitializeEdges(3);
		  Cd->AddEdge(Ei);
		  Cd->AddEdge(Em);
		  Cd->AddEdge(E2);

		  /*
		   * Inizializza e aggiunge i 3 lati della cella di sinistra
		   */
		  Cs->InitializeEdges(3);
		  Cs->AddEdge(Ee);
		  Cs->AddEdge(E1);
		  Cs->AddEdge(Em);

		  /*
		   * Aggiornamento delle celle vicine per i lati Ei ed Ee
		   */
		  Ei->AddCell(Cd);
		  Ee->AddCell(Cs);

		  /*
		   * Gestione del vicinato del nuovo lato creato (Em divide in due il triangolo)
		   * 0:dx, 1:sx
		   */
		  Em->InitializeCells(2); // Genera e aggiunge 2 celle vicine al nuovo Edge
		  Em->AddCell(Cd);
		  Em->AddCell(Cs);


#ifdef DEBUG_NEIGHBORHOOD
		cout<<"\t\tCella figlia destra: Ei-Em-E2\n\t\tCella figlia sinistra: Ee-E1-Em"<<endl;
#endif

		  #ifdef DEBUG_REFINER
		  	  cout<<"\tCreate le figlie: "<<Cd->Id()<<" e "<<Cs->Id()<<endl;
			  cout<<"\t\tCoordinate di "<<Cd->Id()<<": P0 ("<<Cd->Point(0)->X()<<", "<<Cd->Point(0)->Y()<<") - P1 ("<<Cd->Point(1)->X()<<", "<<Cd->Point(1)->Y()<<") - P2 ("<<Cd->Point(2)->X()<<", "<<Cd->Point(2)->Y()<<")"<<endl;
			  cout<<"\t\t\tLati di "<<Cd->Id()<<": "<<Cd->Edge(0)->Id()<<", "<<Cd->Edge(1)->Id()<<" e "<<Cd->Edge(2)->Id()<<endl;

			  /*for(int i=0;i<3;i++){
				  string right;
				  string left;
				  if(Cd->Edge(i)->HasRightCell()){
					  right=to_string(Cd->Edge(i)->RightCell()->Id());
				  } else {
				  	  right="void";
				  }
				  if(Cd->Edge(i)->HasLeftCell()){
				  	  left=to_string(Cd->Edge(i)->LeftCell()->Id());
				  } else {
				  	  left="void";
				  }
				  cout<<"\t\t\t\t\tVicini di "<<Cd->Edge(i)->Id()<<": "<<right<<" e "<<left<<endl;
			  }*/
			  cout<<"\t\tCoordinate di "<<Cs->Id()<<": P0 ("<<Cs->Point(0)->X()<<", "<<Cs->Point(0)->Y()<<") - P1 ("<<Cs->Point(1)->X()<<", "<<Cs->Point(1)->Y()<<") - P2 ("<<Cs->Point(2)->X()<<", "<<Cs->Point(2)->Y()<<")"<<endl;
			  cout<<"\t\t\tLati di "<<Cs->Id()<<": "<<Cs->Edge(0)->Id()<<", "<<Cs->Edge(1)->Id()<<" e "<<Cs->Edge(2)->Id()<<endl;
			  /*for(int i=0;i<3;i++){
				  string right_1;
				  string left_1;
				  if(Cs->Edge(i)->HasRightCell()){
					  right_1=to_string(Cs->Edge(i)->RightCell()->Id());
				  } else {
				  	  right_1="void";
				  }
				  if(Cs->Edge(i)->HasLeftCell()){
					  left_1=to_string(Cs->Edge(i)->LeftCell()->Id());
				  } else {
				  	  left_1="void";
				  }
				  cout<<"\t\t\t\t\tVicini di "<<Cs->Edge(i)->Id()<<": "<<right_1<<" e "<<left_1<<endl;
			  }*/

			  #endif

#ifdef BREAKPOINTS
		  	cout<<"RT BP7"<<endl;
#endif
		  /*
		   * Verifica qual è la cella padre da rimuovere per aggiornare la vicinanza sui due lati E2 ed E1
		   * 0:dx, 1:sx
		   */
		  if(E2->HasRightCell() && E2->RightCell()->Id() == cell->Id() ){
			  E2->InsertCell(Cd,0);
		  }
		  if(E2->HasLeftCell() && E2->LeftCell()->Id() == cell->Id() ){
			  E2->InsertCell(Cd,1);
		  }
		  if(E1->HasRightCell() && E1->RightCell()->Id() == cell->Id() ){
			  E1->InsertCell(Cs,0);
		  }
		  if(E1->HasLeftCell() && E1->LeftCell()->Id() == cell->Id() ){
			  E1->InsertCell(Cs,1);
		  }

#ifdef BREAKPOINTS
		  cout<<"RT BP8"<<endl;
#endif
		  /*
		   * Se uno dei due lati è inattivo (E1 o E2), aggiungo la figlia corrispondente ai triangoli da tagliare, poichè
		   * significa che tale cella è da tagliare per rendere uniforme la mesh
		   */
		  if(!meshPointer->Edge(Cs->Edge(0)->Id())->IsActive() || !meshPointer->Edge(Cs->Edge(1)->Id())->IsActive() || !meshPointer->Edge(Cs->Edge(2)->Id())->IsActive()){
			  AddCellId(Cs->Id());
			 // cout<<"\tAggiungo la cella "<<Cs->Id()<<" (cella figlia sx)"<<endl;
		  }
		  if(!meshPointer->Edge(Cd->Edge(0)->Id())->IsActive() || !meshPointer->Edge(Cd->Edge(1)->Id())->IsActive() || !meshPointer->Edge(Cd->Edge(2)->Id())->IsActive()){
			  AddCellId(Cd->Id());
			 // cout<<"\tAggiungo la cella "<<Cd->Id()<<" (cella figlia dx)"<<endl;
		  }

#ifdef BREAKPOINTS
		  cout<<"RT BP9"<<endl;
#endif
		  return Output::Success;
	  }


	  /*
	   * Funzione per il taglio di un triangolo con tutti e 3 i lati inattivi -> Viene tagliato in 4 triangoli
	   */
	  Output::ExitCodes TriangleRefiner::RefineUniformly( const unsigned int& value ){
		  Rotate(value);
		  GenericCell* cell=meshPointer->Cell(value);
		  vector <GenericEdge*> edges; // lati del triangolo da tagliare
		  vector <GenericPoint*> mediumPoints; // punti medi dei tre lati del triangolo da tagliare
		  GenericEdge *E0,*E1,*E2; // lati che devono essere creati per tagliare il triangolo
		  GenericPoint *C0P0, *C0P1, *C1P0, *C1P1; // Child0Point0...
		  GenericEdge *E0C0, *E0C1, *E1C0, *E1C1, *E2C0, *E2C1;// Edge0

#ifdef BREAKPOINTS
		  cout<<"RU BP1"<<endl;
#endif
		  /*
		   * Ricerca del punto medio per ogni lato: si ricerca il punto in comune tra i due lati figli
		   */
		  for(int i=0; i<3;i++){
			  edges.push_back(meshPointer->Edge(cell->Edge(i)->Id()));
			  C0P0=meshPointer->Point(( (GenericEdge*)(edges[i]->Child(0)) ) ->Point(0)->Id());
			  C0P1=meshPointer->Point(( (GenericEdge*)(edges[i]->Child(0)) ) ->Point(1)->Id());
			  C1P0=meshPointer->Point(( (GenericEdge*)(edges[i]->Child(1)) ) ->Point(0)->Id());
			  C1P1=meshPointer->Point(( (GenericEdge*)(edges[i]->Child(1)) ) ->Point(1)->Id());
			  if(C0P0==C1P0 || C0P0==C1P1){
				  mediumPoints.push_back(C0P0);
			  } else if (C0P1==C1P0 || C0P1==C1P1){
				  mediumPoints.push_back(C0P1);
			  }
		  }

		  /*
		   * Assegnazione dei figli dei lati già esistenti ai puntatori
		   */
		  E0C0=meshPointer->Edge(edges[0]->Child(0)->Id());
		  E0C1=meshPointer->Edge(edges[0]->Child(1)->Id());
		  E1C0=meshPointer->Edge(edges[1]->Child(0)->Id());
		  E1C1=meshPointer->Edge(edges[1]->Child(1)->Id());
		  E2C0=meshPointer->Edge(edges[2]->Child(0)->Id());
		  E2C1=meshPointer->Edge(edges[2]->Child(1)->Id());

		  cell->InitializeChilds(4);

#ifdef BREAKPOINTS
		  cout<<"RU BP2"<<endl;
#endif
		  /*
		   * Creazione e aggiornamento vicinanza per il triangolo T0 (corrispondente
		   * al punto P0 del lato E0).
		   */
		  GenericCell* T0=meshPointer->CreateCell();
		  meshPointer->AddCell(T0);
		  T0->SetFather(cell);
		  cell->AddChild(T0);
		  T0->InitializePoints(3);
		  T0->InitializeEdges(3);
		  T0->AddPoint(edges[0]->Point(0));
		  T0->AddPoint(mediumPoints[0]);
		  T0->AddPoint(mediumPoints[2]);
		  T0->AddEdge(E0C1);
		  E2=meshPointer->CreateEdge();
		  meshPointer->AddEdge(E2);
		  E2->AddPoint(mediumPoints[2]);
		  E2->AddPoint(mediumPoints[0]);
		  T0->AddEdge(E2);
		  T0->AddEdge(E2C0);

		  /*
		   * Creazione e aggiornamento vicinanza per il triangolo T1 (corrispondente
		   * al punto P1 del lato E0).
		   */
		  GenericCell* T1=meshPointer->CreateCell();
		  meshPointer->AddCell(T1);
		  T1->SetFather(cell);
		  cell->AddChild(T1);
		  T1->InitializePoints(3);
		  T1->InitializeEdges(3);
		  T1->AddPoint(mediumPoints[0]);
		  T1->AddPoint(edges[0]->Point(1));
		  T1->AddPoint(mediumPoints[1]);
		  T1->AddEdge(E0C0);
		  E0=meshPointer->CreateEdge();
		  meshPointer->AddEdge(E0);
		  E0->AddPoint(mediumPoints[0]);
		  E0->AddPoint(mediumPoints[1]);
		  T0->AddEdge(E0);
		  T0->AddEdge(E1C1);

		  /*
		   * Creazione e aggiornamento vicinanza per il triangolo T2 (corrispondente
		   * al punto P1 del lato E1).
		   */
		  GenericCell* T2=meshPointer->CreateCell();
		  meshPointer->AddCell(T2);
		  T2->SetFather(cell);
		  cell->AddChild(T2);
		  T2->InitializePoints(3);
		  T2->InitializeEdges(3);
		  T2->AddPoint(mediumPoints[1]);
		  T2->AddPoint(edges[1]->Point(1));
		  T2->AddPoint(mediumPoints[2]);
		  T2->AddEdge(E1C0);
		  E1=meshPointer->CreateEdge();
		  meshPointer->AddEdge(E1);
		  E1->AddPoint(mediumPoints[1]);
		  E1->AddPoint(mediumPoints[2]);
		  T1->AddEdge(E1);
		  T1->AddEdge(E2C1);

		  /*
		   * Creazione e aggiornamento vicinanza per il triangolo T3 (corrispondente
		   * ai 3 punti medi dei lati del triangolo da tagliare).
		   */
		  GenericCell* T3=meshPointer->CreateCell();
		  meshPointer->AddCell(T3);
		  T3->SetFather(cell);
		  cell->AddChild(T3);
		  T3->InitializePoints(3);
		  T3->InitializeEdges(3);
		  T3->AddPoint(mediumPoints[0]);
		  T3->AddPoint(mediumPoints[1]);
		  T3->AddPoint(mediumPoints[2]);
		  T3->AddEdge(E0);
		  T3->AddEdge(E1);
		  T3->AddEdge(E2);

		  cell->SetState(false);

#ifdef BREAKPOINTS
		  cout<<"RU BP3"<<endl;
#endif
		  /*
		   * Aggiornamento vicinanza per i lati esterni
		   */
		  E0C1->InitializeCells(2);
		  E0C1->AddCell(T0);

		  E0C0->InitializeCells(2);
		  E0C0->AddCell(T1);

		  E1C1->InitializeCells(2);
		  E1C1->AddCell(T1);

		  E1C0->InitializeCells(2);
		  E1C0->AddCell(T2);

		  E2C1->InitializeCells(2);
		  E2C1->AddCell(T2);

		  E2C0->InitializeCells(2);
		  E2C0->AddCell(T0);

		  /*
		   * Aggiornamento vicinanza per i lati interni
		   */
		  E0->InitializeCells(2);
		  E0->AddCell(T1);
		  E0->AddCell(T3);

		  E1->InitializeCells(2);
		  E1->AddCell(T3);
		  E1->AddCell(T2);

		  E2->InitializeCells(2);
		  E2->AddCell(T0);
		  E2->AddCell(T3);

#ifdef BREAKPOINTS
		  cout<<"RU BP4"<<endl;
#endif
		  return Output::Success;
	  }

	  /*
	   * Aggiornare la vicinanza anche per i punti e per le celle (point-edge, point-cell, cell-cell)
	   * Da mettere nella relazione
	   */
	  void TriangleRefiner::UpdateNeighbourhood(){
		  GenericEdge* lato, *lati[2];
		  GenericPoint* punto, *p;
		  GenericCell *cella;
		  unsigned int i, j, k, h, t, n;
		  bool trovato, trovato1, trovato2;

		  for(i=0; i<meshPointer->NumberOfCells(); i++){
			  cella = meshPointer->Cell(i);
			  if(cella->IsActive()){  /* Completa la vicinanza solo se la cella è attiva */
				  	for(j=0; j<3; j++){ /* Per i 3 punti ed i 3 lati della cella */
					  lato = meshPointer->Edge(cella->Edge(j)->Id());
					  punto = meshPointer->Point(cella->Point(j)->Id());

					  /*
					   * Vicinanza cell->cell: Aggiunge la cella vicina trovata ai vicini della cella che si sta elaborando
					   *
					   */
					  trovato=false;
					  /* Se la cella destra del lato j-esimo non è la cella attuale, verifico se devo aggiungerla*/
					  if(lato->HasRightCell() && lato->RightCell()->Id() != cella->Id()){
						  /* Scorre il vettore delle celle vicine, se la cella dx del lato che si sta considerando non è ancora presente la aggiunge */
						  for(k=0; k < cella->NumberOfCells(); k++){
							  /* Se trova una cella inattiva la rimuove dalle celle vicine */
							  if(cella->Cell(k)!=NULL &&  !cella->Cell(k)->IsActive()){
								  cella->EraseCell(k);
								  //cout<<"Rimuovo cella inattiva (cella)"<<endl;
							  }
							  if(cella->Cell(k)!=NULL &&  cella->Cell(k)->Id() == lato->RightCell()->Id()){
								  trovato=true;
							  }
						  }
						  if(!trovato && lato->RightCell()->IsActive()){
							  cella->AddCell(lato->RightCell());
						  }

					  } else if(lato->HasLeftCell() && lato->LeftCell()->Id() != cella->Id()){
						  for(k=0; k < cella->NumberOfCells(); k++){
							  if(cella->Cell(k)!=NULL &&  !cella->Cell(k)->IsActive()){
								  cella->EraseCell(k);
								//  cout<<"Rimuovo cella inattiva (cella)"<<endl;
							  }
							  if(cella->Cell(k)!=NULL &&  cella->Cell(k)->Id() == lato->LeftCell()->Id()){
								  trovato=true;
							  }
						  }
						  if(!trovato && lato->LeftCell()->IsActive()){
							  cella->AddCell(lato->LeftCell());
						  }
					  }


					  /*
					   * Vicinanza point->cell: Aggiunge la cella che si sta elaborando come vicina del punto
					   *
					   */
					  trovato=false;
					  for(k=0; k<punto->NumberOfCells(); k++){
						  /* Se trova una cella inattiva tra le vicinanze del punto la rimuove*/
						  if(punto->Cell(k)!=NULL &&  !punto->Cell(k)->IsActive()){
							  punto->EraseCell(k);
							  //cout<<"Rimuovo cella inattiva (punto)"<<endl;
						  }
						  /* Se il punto che si sta considerando non contiene ancora l'informazione di vicinanza con la cella attuale, la aggiunge */
						  if(punto->Cell(k)!=NULL && punto->Cell(k)->Id() == cella->Id()){
							  trovato=true;
						  }
					  }
					  if(!trovato){
						  punto->AddCell(cella);
					  }


					  /*
					   * Vicinanza point->edge: Aggiunge il lato ai vicini dei suoi due punti
					   *
					   */
					  h=0;
					  while(h<2){
						  p=meshPointer->Point(lato->Point(h)->Id());
						  trovato = false;
						  for(k=0;k < p->NumberOfEdges(); k++){
							  /* Se trova un lato inattivo tra le vicinanze del punto lo rimuove */
						  	  if(p->Edge(k)!=NULL &&  !p->Edge(k)->IsActive()){
						  		  p->EraseEdge(k);
						  		  //cout<<"Rimuovo lato inattivo (punto)"<<endl;
						  	  }
						  	  if( p->Edge(k)->Id() == lato->Id() ){
								  trovato=true;
							  }
						  }
						  if(!trovato && lato->IsActive()){
							  p->AddEdge(lato);
						  }
						  h++;
					  }
				  } // for(...j<3...)
			  } // if(cella->IsActive())
		  } // for(..i<meshPointer->NumberOfCells()..)
		  return;
	  }

	  /*
	   * Funzione di debug per la stampa delle vicinanze dei punti
	   * Non serve nella relazione
	   */
	  void TriangleRefiner::PrintNeigh(){
		  for(unsigned int i=0; i<meshPointer->NumberOfPoints(); i++){
			  if(meshPointer->Point(i)->IsActive()){
				  cout<<"\nPunto "<<i<<" ("<<meshPointer->Point(i)->X()<<", "<<meshPointer->Point(i)->Y()<<")\n\t\tcelle vicine: ";
				  for(unsigned int j=0; j<meshPointer->Point(i)->NumberOfCells(); j++){
						 cout<<meshPointer->Point(i)->Cell(j)->Id()<<", ";
				  }
				  cout<<endl;
				  cout<<"\t\tlati vicini: ";
				  for(unsigned int j=0; j<meshPointer->Point(i)->NumberOfEdges(); j++){
						 cout<<meshPointer->Point(i)->Edge(j)->Id()<<", ";
				  }
				  cout<<endl;
			  }
		  }
	  }


	  /*
	   * Initialize* vector::reserve -> Requests that the vector capacity be at least enough to contain n elements. If n is greater than the current vector capacity, the function causes the container to reallocate its storage increasing its capacity to n (or greater).
	   * 	In all other cases, the function call does not cause a reallocation and the vector capacity is not affected. This function has no effect on the vector size and cannot alter its elements.
	   */

	  /*
	   * Allocate* vector::resize -> Resizes the container so that it contains n elements. If n is smaller than the current container size, the content is reduced to its first n elements, removing those beyond (and destroying them).
	   * 	If n is greater than the current container size, the content is expanded by inserting at the end as many elements as needed to reach a size of n. If val is specified, the new elements are initialized as copies of val, otherwise, they are value-initialized.
	   * 	If n is also greater than the current container capacity, an automatic reallocation of the allocated storage space takes place.
	   */
}
