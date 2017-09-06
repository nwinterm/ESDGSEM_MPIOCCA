
#include "MeshPartitioning.h"


using namespace std;
MeshPartitioning::MeshPartitioning( const int procs, const int N)
{
NumProcessors=procs;
    ngl = N;
    ngl2 = N*N;
    CommTags.resize(NumProcessors,NumProcessors);
    int counter=0;
    for (int i=1;i<=NumProcessors;i++){
        for (int j=1;j<=NumProcessors;j++){
            CommTags(i,j) = counter;
            counter++;
        }
    }
}

MeshPartitioning::~MeshPartitioning()
{
    //dtor
}



void MeshPartitioning::DivideMesh(const Mesh GlobalMesh,const MPI_setup MPI)
{



   global_NumElements=GlobalMesh.m_num_elements;

   	 for (int i=1;i<NumProcessors;i++){
        MPI_Isend (&global_NumElements,1,MPI_INT,i,i,MPI_COMM_WORLD,&MPI.reqs[i]);
	 }
    MPI_Waitall(NumProcessors,MPI.reqs, MPI.stats);

   global_NumEdges=GlobalMesh.m_num_edges;

    ElementsPerProc.resize(NumProcessors);
    EdgesPerProc.resize(NumProcessors);


      cout <<"using " << NumProcessors <<" processors for a problem with "<< global_NumElements <<" Elements.\n";
      for (int i=1;i<=NumProcessors;i++){
        ElementsPerProc(i) = global_NumElements / NumProcessors;
      }

      int RestElements = global_NumElements % NumProcessors;

	 for (int i=1;i<=RestElements;i++){
			ElementsPerProc(i)=ElementsPerProc(i)+1;
	 }


	 NumElements=ElementsPerProc(1); // set number of elements for host
	 //SEND number of elements to other tasks
	 for (int i=1;i<NumProcessors;i++){
            cout << "Elements in processor " << i+1 <<" : " << ElementsPerProc(i+1) <<"\n";
            MPI_Isend (&ElementsPerProc(i+1),1,MPI_INT,i,i,MPI_COMM_WORLD,&MPI.reqs[i]);

	 }
    MPI_Waitall(NumProcessors,MPI.reqs, MPI.stats);





    ElementLocalToGlobal.resize(ElementsPerProc(1),NumProcessors);  //Link LOCAL Elements for EACH Processors to GLOBAL ELEMENT ID
    ElementGlobalToLocal.resize(2,global_NumElements); //Link each GLOBAL element to LOCALELEMENT + PROCESSOR

    int eleID_global=1;
    for (int iproc=1;iproc<=NumProcessors;iproc++){
        for (int ie=1;ie<=ElementsPerProc(iproc);ie++){

			ElementLocalToGlobal(ie,iproc) = eleID_global;
			ElementGlobalToLocal(1,eleID_global) = ie;
			ElementGlobalToLocal(2,eleID_global) = iproc;
            eleID_global++;

        }
    }

    MyElementLocalToGlobal.resize(NumElements);
    for( int ie = 1; ie<=ElementsPerProc(1);ie++){
        MyElementLocalToGlobal(ie) = ElementLocalToGlobal(ie,1);
    }

    for (int iproc=1;iproc<NumProcessors;iproc++){
            MPI_Isend (&ElementLocalToGlobal(1,iproc+1),ElementsPerProc(iproc+1),MPI_INT,iproc,iproc,MPI_COMM_WORLD,&MPI.reqs[iproc]);
//        MyElementLocalToGlobal(ie) = ElementLocalToGlobal(ie,1)
    }

    MPI_Waitall(NumProcessors,MPI.reqs, MPI.stats);


//cout << "ELEMENT SPLITTING FINISHED! \n" ;


// ///////////////////////////////////////////////////////////////////////////  EDGES ///////////////////////////////////////////////

//    for (int i=1;i<=global_NumElements;i++){
//     cout << "global ID: " << i << "\n";
//     cout << " Side 1: " << GlobalMesh.ElementToEdge(1,i) << " Side 2: " << GlobalMesh.ElementToEdge(2,i) <<" Side 3: " << GlobalMesh.ElementToEdge(3,i) <<" Side 4: " << GlobalMesh.ElementToEdge(4,i) <<"\n";
//    }


int * EdgeArray_proc = (int*) calloc(global_NumEdges,sizeof(int));

//int EdgeArray_proc[global_NumEdges];
int LocatedAt;
int maxEdgesLocal=0;
    EdgeLocalToGlobal.resize(global_NumEdges,NumProcessors);
    EdgeGlobalToLocal.resize(4,global_NumEdges);

    for (int i=1;i<=global_NumEdges;i++){
            EdgeArray_proc[i-1]=-1;
            EdgeGlobalToLocal(1,i)=0;
            EdgeGlobalToLocal(2,i)=0;
            EdgeGlobalToLocal(3,i)=0;
            EdgeGlobalToLocal(4,i)=0;
        for (int iproc=1;iproc<=NumProcessors;iproc++){
            EdgeLocalToGlobal(i,iproc)=-1;
        }
    }
    for (int iproc=1;iproc<=NumProcessors;iproc++){
        EdgesPerProc(iproc) = 0;
    }

    for (int iproc=1;iproc<=NumProcessors;iproc++){

        // use an array to store all global edges that belong to this processor
        for (int i=1;i<=global_NumEdges;i++){
            EdgeArray_proc[i-1] = -1; //EdgeLocalToGlobal(i,iproc);
        }

        for (int ie=1;ie<=ElementsPerProc(iproc);ie++){
            eleID_global = ElementLocalToGlobal(ie,iproc);

            for (int is=1;is<=4;is++){
                //get global side ID
                int sideID_global=GlobalMesh.ElementToEdge(is,eleID_global);



//                cout << "eleID_global: " <<eleID_global <<"\n";
//                cout << "sideID_global: " <<sideID_global <<"\n";

                // CHECK IF WE ALREADY FOUND THIS EDGE ON THIS PROCESSOR
                isInArray(sideID_global,EdgeArray_proc,EdgesPerProc(iproc),&LocatedAt);
//                cout << "LocatedAt: " <<LocatedAt <<"\n";
//                if (sideID_global == 3316) {
//                    cout << "EDGE: " << sideID_global <<" was found at "<< LocatedAt <<"\n";
//
//                }
//                cout << "EdgeInfo(3,sideID_global)" << GlobalMesh.EdgeInfo(3,sideID_global)+1 <<"\n";
                if (LocatedAt == 0) {
                    //PRINT*, 'SIDE NOT FOUND: ', sideID_global
                    EdgesPerProc(iproc) = EdgesPerProc(iproc)+1;
//                    cout << " EDGE ADDED : " << sideID_global << "Element: " << eleID_global << " side " << is <<"\n" ;

                    EdgeLocalToGlobal(EdgesPerProc(iproc),iproc) = sideID_global;
                    EdgeArray_proc[EdgesPerProc(iproc)-1] = sideID_global; //EdgeLocalToGlobal(i,iproc);
                    //check if this element is LEFT or RIGHT element to the edge
                    if (GlobalMesh.EdgeInfo(3,sideID_global)==eleID_global-1){
                          EdgeGlobalToLocal(1,sideID_global) = EdgesPerProc(iproc);
                          EdgeGlobalToLocal(3,sideID_global) = iproc;
                    }else{

                          EdgeGlobalToLocal(2,sideID_global) = EdgesPerProc(iproc);
                          EdgeGlobalToLocal(4,sideID_global) = iproc;
                    }
                }else{

                    if (GlobalMesh.EdgeInfo(3,sideID_global)==eleID_global-1){
//                        cout << " WE ACTUALLY GET IN WEIRD CASE WHERE RIGHT ELE WAS FOUND FIRST\n";
                          EdgeGlobalToLocal(1,sideID_global) = EdgeGlobalToLocal(2,sideID_global);
                          EdgeGlobalToLocal(3,sideID_global) = EdgeGlobalToLocal(4,sideID_global);
                    }else{
//                        if (sideID_global>3280){cout <<"STIMMT WAS NICHT!";}
                          EdgeGlobalToLocal(2,sideID_global) = EdgeGlobalToLocal(1,sideID_global);
                          EdgeGlobalToLocal(4,sideID_global) = EdgeGlobalToLocal(3,sideID_global);
                    }
                }

            }

        }
        maxEdgesLocal = max(maxEdgesLocal,EdgesPerProc(iproc))	;
    }


//    for (int i=global_NumEdges-35;i<=global_NumEdges;i++){
//     cout << "global ID: " << i << "\n";
//     cout << "Proc Left: "<< EdgeGlobalToLocal(3,i) << " Proc Right: " << EdgeGlobalToLocal(4,i) <<"\n";
//     cout << " Local Side Left: " << EdgeGlobalToLocal(1,i) <<  " LocalSide Right: " << EdgeGlobalToLocal(2,i) <<"\n";
//
//
//    }

//cout << "EDGE SPLITTING FINISHED! Max Edges Local are "<<  maxEdgesLocal  <<"\n" ;



int * globalEdgeInfo = (int*) calloc(10*maxEdgesLocal*NumProcessors,sizeof(int));

int localEdgeCounter[NumProcessors];

for (int iproc=0;iproc<NumProcessors;iproc++){
    localEdgeCounter[iproc] = 0;
}



for (int is=1;is<=global_NumEdges;is++){
  int edgeIDlocal_LEFT  = EdgeGlobalToLocal(1,is);		// 1=LOCAL SIDE ID PROC LEFT,  2=LOCAL SIDE ID PROC RIGHT,  3=LEFT PROCESSOR, 4=RIGHT PROCESSOR
  int edgeIDlocal_RIGHT = EdgeGlobalToLocal(2,is);
  int proc_LEFT	    = EdgeGlobalToLocal(3,is);
  int proc_RIGHT	    = EdgeGlobalToLocal(4,is);
  localEdgeCounter[proc_LEFT-1]++;


//  cout << "global edge id: " << is <<"\n";
//  cout << "edgeIDlocal_LEFT: " << edgeIDlocal_LEFT <<"\n";
//  cout << "edgeIDlocal_RIGHT: " << edgeIDlocal_RIGHT <<"\n";
//  cout << "proc_LEFT: " << proc_LEFT <<"\n";
//  cout << "proc_RIGHT: " << proc_RIGHT <<"\n";

//CHECK IF THIS EDGE CROSSES PROCESSOR LINES
	if (proc_LEFT != proc_RIGHT) {
		if (proc_RIGHT == 0 ) {
            int id1 = 10*maxEdgesLocal*(proc_LEFT-1)+ 10*(edgeIDlocal_LEFT-1);
            globalEdgeInfo[id1+0] = GlobalMesh.EdgeInfo(1,is);//				!these are not even used after the inital mesh generation
			globalEdgeInfo[id1+1] = GlobalMesh.EdgeInfo(2,is);//			! ""
			globalEdgeInfo[id1+2] = ElementGlobalToLocal(1,GlobalMesh.EdgeInfo(3,is)+1)-1	;//this is changed to be the local element on left processor
			globalEdgeInfo[id1+3] = -1; //neighbour element is zero since its a boundary
			globalEdgeInfo[id1+4] = GlobalMesh.EdgeInfo(5,is)	;//local side left element <- fine as is
			globalEdgeInfo[id1+5] = GlobalMesh.EdgeInfo(6,is);	//local side right element <- fine as is
			globalEdgeInfo[id1+6] = GlobalMesh.EdgeInfo(7,is);
			globalEdgeInfo[id1+7] = proc_LEFT -1;
			globalEdgeInfo[id1+8] = proc_RIGHT -1;
			globalEdgeInfo[id1+9] = is -1	;	//global side id
		}else{
			localEdgeCounter[proc_RIGHT-1]++;

            int id1 = 10*maxEdgesLocal*(proc_LEFT-1)+ 10*(edgeIDlocal_LEFT-1);
            int id2 = 10*maxEdgesLocal*(proc_RIGHT-1)+ 10*(edgeIDlocal_RIGHT-1);

			globalEdgeInfo[id1+0] = GlobalMesh.EdgeInfo(1,is);				//these are not even used after the inital mesh generation
			globalEdgeInfo[id1+1] = GlobalMesh.EdgeInfo(2,is)		;		//! ""
			globalEdgeInfo[id1+2] = ElementGlobalToLocal(1,GlobalMesh.EdgeInfo(3,is)+1)-1	;//!this is changed to be the local element on left processor
			globalEdgeInfo[id1+3] = GlobalMesh.EdgeInfo(4,is) 	;//!this is changed to be the GLOBAL ELEMENT ID of the right element!!!!
			globalEdgeInfo[id1+4] = GlobalMesh.EdgeInfo(5,is)	;//!local side left element <- fine as is
			globalEdgeInfo[id1+5] = GlobalMesh.EdgeInfo(6,is)	;//!local side right element <- fine as is
			globalEdgeInfo[id1+6] = GlobalMesh.EdgeInfo(7,is);
			globalEdgeInfo[id1+7] = proc_LEFT -1;
			globalEdgeInfo[id1+8] = proc_RIGHT -1;
			globalEdgeInfo[id1+9] = is -1 	;	//!global side id

			globalEdgeInfo[id2+0] = GlobalMesh.EdgeInfo(1,is)	;	//		!these are not even used after the inital mesh generation
			globalEdgeInfo[id2+1]  = GlobalMesh.EdgeInfo(2,is);			//	! ""
			globalEdgeInfo[id2+2]  = GlobalMesh.EdgeInfo(3,is),	//!this is changed to be the GLOBAL ELEMENT ID on left processor
			globalEdgeInfo[id2+3]  = ElementGlobalToLocal(1,GlobalMesh.EdgeInfo(4,is)+1)-1;//	!this is changed to be the local element on right processor
			globalEdgeInfo[id2+4]  = GlobalMesh.EdgeInfo(5,is);	//!local side left element <- fine as is
			globalEdgeInfo[id2+5]  = GlobalMesh.EdgeInfo(6,is);//!local side right element <- fine as is
			globalEdgeInfo[id2+6]  = GlobalMesh.EdgeInfo(7,is);
			globalEdgeInfo[id2+7]  = proc_LEFT -1;
			globalEdgeInfo[id2+8]  = proc_RIGHT -1;
			globalEdgeInfo[id2+9]  = is -1	;	//!global side id
        }
	}else{

            int id1 = 10*maxEdgesLocal*(proc_LEFT-1)+ 10*(edgeIDlocal_LEFT-1);
			globalEdgeInfo[id1+0]= GlobalMesh.EdgeInfo(1,is)	;		//	!these are not even used after the inital mesh generation
			globalEdgeInfo[id1+1] = GlobalMesh.EdgeInfo(2,is)	;		//	! ""
			globalEdgeInfo[id1+2] = ElementGlobalToLocal(1,GlobalMesh.EdgeInfo(3,is)+1)-1	; //!this is changed to be the local element on left processor
			globalEdgeInfo[id1+3] = ElementGlobalToLocal(1,GlobalMesh.EdgeInfo(4,is)+1)-1 ;	//!this is changed to be the local element on right processor
			globalEdgeInfo[id1+4] = GlobalMesh.EdgeInfo(5,is);	//!local side left element <- fine as is
			globalEdgeInfo[id1+5] = GlobalMesh.EdgeInfo(6,is)	;//!local side right element <- fine as is
			globalEdgeInfo[id1+6] = GlobalMesh.EdgeInfo(7,is);
			globalEdgeInfo[id1+7]= proc_LEFT -1;
			globalEdgeInfo[id1+8]= proc_LEFT -1;
			globalEdgeInfo[id1+9] = is -1	;	//!global side id

    }
}



NumEdges=EdgesPerProc(1);// set edge number for host

//cout << "My NumEdges is " << NumEdges <<"\n";
cout << "Edges in processor " << 1 <<" : " << EdgesPerProc(1) <<"\n";
 for (int i=1;i<NumProcessors;i++){
        cout << "Edges in processor " << i+1 <<" : " << EdgesPerProc(i+1) <<"\n";
        MPI_Isend (&EdgesPerProc(i+1),1,MPI_INT,i,i,MPI_COMM_WORLD,&MPI.reqs[i]);
 }
MPI_Waitall(NumProcessors,MPI.reqs, MPI.stats);
  for (int i=1;i<NumProcessors;i++){
        MPI_Isend (&global_NumEdges,1,MPI_INT,i,i,MPI_COMM_WORLD,&MPI.reqs[i]);
 }
MPI_Waitall(NumProcessors,MPI.reqs, MPI.stats);

MyEdgesLocalToGlobal.resize(NumEdges);
MyEdgeInfo.resize(10,NumEdges);

for (int i=1;i<=EdgesPerProc(1);i++){
   int id1 = i-1;

    MyEdgesLocalToGlobal(i) = EdgeLocalToGlobal(i,1);

    int id =  10*(i-1);
    for (int j=1;j<=10;j++){
        MyEdgeInfo(j,i)= globalEdgeInfo[id+j-1];
    }

}


    for (int iproc=1;iproc<NumProcessors;iproc++){
        MPI_Isend (&EdgeLocalToGlobal(1,iproc+1),EdgesPerProc(iproc+1),MPI_INT,iproc,iproc,MPI_COMM_WORLD,&MPI.reqs[iproc]);
    }
    MPI_Waitall(NumProcessors,MPI.reqs, MPI.stats);

    for (int iproc=1;iproc<NumProcessors;iproc++){
        imatrix globalEdgeInfo_TMP(10,EdgesPerProc(iproc+1));
        for (int i=1;i<=EdgesPerProc(iproc+1);i++){
            int id = 10*maxEdgesLocal*iproc+ 10*(i-1);
            for (int j=1;j<=10;j++){
                globalEdgeInfo_TMP(j,i)= globalEdgeInfo[id+j-1];
            }

        }
        MPI_Isend (&globalEdgeInfo_TMP(1,1),10*EdgesPerProc(iproc+1),MPI_INT,iproc,iproc,MPI_COMM_WORLD,&MPI.reqs[iproc]);
        MPI_Wait(&MPI.reqs[iproc], MPI.stats);
    }




x_global.resize(ngl*ngl*NumElements,1);
SplitRealValuesBetweenProcs( MPI,GlobalMesh.x_global , x_global);


y_global.resize(ngl*ngl*NumElements,1);
SplitRealValuesBetweenProcs( MPI,GlobalMesh.y_global , y_global);

yEta_global.resize(ngl*ngl*NumElements,1);
SplitRealValuesBetweenProcs( MPI,GlobalMesh.yEta_global , yEta_global);
xEta_global.resize(ngl*ngl*NumElements,1);
SplitRealValuesBetweenProcs( MPI,GlobalMesh.xEta_global , xEta_global);
yXi_global.resize(ngl*ngl*NumElements,1);
SplitRealValuesBetweenProcs( MPI,GlobalMesh.yXi_global , yXi_global);
xXi_global.resize(ngl*ngl*NumElements,1);
SplitRealValuesBetweenProcs( MPI,GlobalMesh.xXi_global , xXi_global);
J_global.resize(ngl*ngl*NumElements,1);
SplitRealValuesBetweenProcs( MPI,GlobalMesh.J_global , J_global);


//    for(int ie=0;ie<NumElements;++ie){
//    cout <<"\nELE: " << ie <<"\n";
//      for(int j=0;j<ngl;++j){
//        for(int i=0;i<ngl;++i){
//            int id = ie*ngl*ngl   +j*ngl+i +1;
//
//            cout << J_global(id,1) << " ";
//
//      }
//       cout  <<"\n";
//    }
//}



//cout  <<" EDGE STUFF   \n";
nx_global.resize(ngl*NumEdges,1);
SplitEdgeRealValuesBetweenProcs( MPI,GlobalMesh.NormalsX , nx_global);

//    for(int ie=0;ie<NumEdges;++ie){
//    cout <<"\n EDGE: " << ie <<"\n";
//      for(int j=0;j<ngl;++j){
//            int id = ie*ngl   +j +1;
//
//            cout << nx_global(id,1) << " ";
//
//
//       cout  <<"\n";
//    }
//}

ny_global.resize(ngl*NumEdges,1);
SplitEdgeRealValuesBetweenProcs( MPI,GlobalMesh.NormalsY , ny_global);
scal_global.resize(ngl*NumEdges,1);
SplitEdgeRealValuesBetweenProcs( MPI,GlobalMesh.Scal , scal_global);





int * SplitElementToEdge = (int*) calloc(4*ElementsPerProc(1)*NumProcessors,sizeof(int));
int * SplitElemEdgeMasterSlave = (int*) calloc(4*ElementsPerProc(1)*NumProcessors,sizeof(int));
//int SplitElementToEdge[4][ElementsPerProc(1)][NumProcessors];
//int SplitElemEdgeMasterSlave[4][ElementsPerProc(1)][NumProcessors];

for (int ie=1;ie<=global_NumElements;ie++){

    int eleID_local=ElementGlobalToLocal(1,ie);
    int iproc= ElementGlobalToLocal(2,ie);

    for (int is=1;is<=4;is++){
        int edgeID = GlobalMesh.ElementToEdge(is,ie);
        int id = 4*(iproc-1)*ElementsPerProc(1) + 4*(eleID_local-1) + (is-1);
        if (GlobalMesh.EdgeInfo(3,edgeID) == ie-1){
            SplitElementToEdge[id] = EdgeGlobalToLocal(1,edgeID);
//            SplitElementToEdge[is-1][eleID_local-1][iproc-1]  = EdgeGlobalToLocal(1,edgeID);
        }else{
            SplitElementToEdge[id] = EdgeGlobalToLocal(2,edgeID);
//            SplitElementToEdge[is-1][eleID_local-1][iproc-1]  = EdgeGlobalToLocal(2,edgeID);
        }
        SplitElemEdgeMasterSlave[id] = GlobalMesh.ElemEdgeMasterSlave(is,ie);

    }
}


ElementToEdge.resize(4,ElementsPerProc(1));
ElemEdgeMasterSlave.resize(4,ElementsPerProc(1));





for (int ie=1;ie<=ElementsPerProc(1);ie++){
    for (int is=1;is<=4;is++){
        int id =  4*(ie-1) + (is-1);
        ElementToEdge(is,ie) = SplitElementToEdge[id];
//        ElementToEdge(is,ie) = SplitElementToEdge[is-1][ie-1][0];
        ElemEdgeMasterSlave(is,ie) = SplitElemEdgeMasterSlave[id];
    }
}

for (int iproc=2;iproc<=NumProcessors;iproc++){
    imatrix EleToEdge_TMP(4,ElementsPerProc(iproc));
    imatrix EleToEdgeMS_TMP(4,ElementsPerProc(iproc));
    for (int ie=1;ie<=ElementsPerProc(iproc);ie++){
        for (int is=1;is<=4;is++){
             int id = 4*(iproc-1)*ElementsPerProc(1) + 4*(ie-1) + (is-1);
            EleToEdge_TMP(is,ie) = SplitElementToEdge[id];
//            EleToEdge_TMP(is,ie) = SplitElementToEdge[is-1][ie-1][iproc-1];
            EleToEdgeMS_TMP(is,ie) = SplitElemEdgeMasterSlave[id];
        }
    }

     MPI_Isend (&EleToEdge_TMP(1,1),4*ElementsPerProc(iproc),MPI_INT,iproc-1,iproc-1,MPI_COMM_WORLD,&MPI.reqs[iproc-1]);
    MPI_Wait(&MPI.reqs[iproc-1], MPI.stats);
     MPI_Isend (&EleToEdgeMS_TMP(1,1),4*ElementsPerProc(iproc),MPI_INT,iproc-1,iproc-1,MPI_COMM_WORLD,&MPI.reqs[iproc-1]);
     MPI_Wait(&MPI.reqs[iproc-1], MPI.stats);
}


















free(EdgeArray_proc);
free(globalEdgeInfo);
free(SplitElementToEdge);
free(SplitElemEdgeMasterSlave);


}

void MeshPartitioning::SplitRealValuesBetweenProcs(const MPI_setup MPI,const fmatrix Values, fmatrix &HostValues){


    for (int ie=1;ie<=ElementsPerProc(1);ie++){
        for (int i=1;i<=ngl;i++){
            for (int j=1;j<=ngl;j++){
                int id = (ElementLocalToGlobal(ie,1)-1)*ngl*ngl   +(j-1)*ngl+i;
                int idLoc = (ie-1)*ngl*ngl   +(j-1)*ngl+i;
                HostValues(idLoc,1) = Values(id,1);
            }
        }
    }





for (int iproc=2;iproc<=NumProcessors;iproc++){
    fmatrix tmpValues(ngl*ngl*ElementsPerProc(iproc));
    for (int ie=1;ie<=ElementsPerProc(iproc);ie++){
        for (int i=1;i<=ngl;i++){
            for (int j=1;j<=ngl;j++){
                int id = (ElementLocalToGlobal(ie,iproc)-1)*ngl*ngl   +(j-1)*ngl+i;
                int idLoc = (ie-1)*ngl*ngl   +(j-1)*ngl+i;
                tmpValues(idLoc,1) = Values(id,1);
            }
        }
    }

        MPI_Isend (&tmpValues(1,1),ngl*ngl*ElementsPerProc(iproc),MPI_DOUBLE,iproc-1,iproc-1,MPI_COMM_WORLD,&MPI.reqs[iproc-1]);
        MPI_Wait(&MPI.reqs[iproc-1], MPI.stats);


}
}

void MeshPartitioning::SplitEdgeRealValuesBetweenProcs(const MPI_setup MPI,const fmatrix Values, fmatrix &HostValues){


    for (int ie=1;ie<=EdgesPerProc(1);ie++){
        for (int i=1;i<=ngl;i++){
                int id = (EdgeLocalToGlobal(ie,1)-1)*ngl   +i;
                int idLoc = (ie-1)*ngl  +i;
                HostValues(idLoc,1) = Values(id,1);
        }
    }


//    for(int ie=0;ie<EdgesPerProc(1);++ie){
//    cout <<"\n EDGE: " << ie <<"\n";
//      for(int j=0;j<ngl;++j){
//            int id = ie*ngl   +j +1;
//
//            cout << HostValues(id,1) << " ";
//
//
//       cout  <<"\n";
//    }
//}



for (int iproc=2;iproc<=NumProcessors;iproc++){
    fmatrix tmpValues(ngl*EdgesPerProc(iproc));
    for (int ie=1;ie<=EdgesPerProc(iproc);ie++){
        for (int i=1;i<=ngl;i++){
                int id = (EdgeLocalToGlobal(ie,iproc)-1)*ngl   +i;
                int idLoc = (ie-1)*ngl  +i;
                tmpValues(idLoc,1) = Values(id,1);
        }
    }

        MPI_Isend (&tmpValues(1,1),ngl*EdgesPerProc(iproc),MPI_DOUBLE,iproc-1,iproc-1,MPI_COMM_WORLD,&MPI.reqs[iproc-1]);
        MPI_Wait(&MPI.reqs[iproc-1], MPI.stats);


}
}




void MeshPartitioning::ReceiveMesh(const MPI_setup MPI){

    MPI_Irecv(&global_NumElements,1,MPI_INT,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);
//    MPI_Wait(&reqs[MPI.rank], stats);
    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);


    MPI_Irecv(&NumElements,1,MPI_INT,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);
//    MPI_Wait(&reqs[MPI.rank], stats);
    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);
//    MPI_Waitall(NumProcessors,MPI.reqs, MPI.stats);

    MyElementLocalToGlobal.resize(NumElements);


    MPI_Irecv(&MyElementLocalToGlobal(1),NumElements,MPI_INT,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);

    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);



    MPI_Irecv(&NumEdges,1,MPI_INT,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);
//    MPI_Wait(&reqs[MPI.rank], stats);
    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);

    MPI_Irecv(&global_NumEdges,1,MPI_INT,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);
//    MPI_Wait(&reqs[MPI.rank], stats);
    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);

//    cout << " IAM RANK: " << MPI.rank << " AND THER ARE " << global_NumEdges << " EDGES GLOBALLY!\n";


    MyEdgeInfo.resize(10,NumEdges);
    MyEdgesLocalToGlobal.resize(NumEdges);
    MPI_Irecv(&MyEdgesLocalToGlobal(1),NumEdges,MPI_INT,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);
    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);

    MPI_Irecv(&MyEdgeInfo(1,1),10*NumEdges,MPI_INT,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);
    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);


    x_global.resize(ngl*ngl*NumElements,1);
    MPI_Irecv(&x_global(1,1),ngl*ngl*NumElements,MPI_DOUBLE,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);
    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);

    y_global.resize(ngl*ngl*NumElements,1);
    MPI_Irecv(&y_global(1,1),ngl*ngl*NumElements,MPI_DOUBLE,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);
    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);

    yEta_global.resize(ngl*ngl*NumElements,1);
    MPI_Irecv(&yEta_global(1,1),ngl*ngl*NumElements,MPI_DOUBLE,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);
    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);

    xEta_global.resize(ngl*ngl*NumElements,1);
    MPI_Irecv(&xEta_global(1,1),ngl*ngl*NumElements,MPI_DOUBLE,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);
    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);

    yXi_global.resize(ngl*ngl*NumElements,1);
    MPI_Irecv(&yXi_global(1,1),ngl*ngl*NumElements,MPI_DOUBLE,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);
    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);

    xXi_global.resize(ngl*ngl*NumElements,1);
    MPI_Irecv(&xXi_global(1,1),ngl*ngl*NumElements,MPI_DOUBLE,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);
    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);

    J_global.resize(ngl*ngl*NumElements,1);
    MPI_Irecv(&J_global(1,1),ngl*ngl*NumElements,MPI_DOUBLE,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);
    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);



nx_global.resize(ngl*NumEdges,1);
MPI_Irecv(&nx_global(1,1),ngl*NumEdges,MPI_DOUBLE,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);
    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);
ny_global.resize(ngl*NumEdges,1);
MPI_Irecv(&ny_global(1,1),ngl*NumEdges,MPI_DOUBLE,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);
    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);
scal_global.resize(ngl*NumEdges,1);
MPI_Irecv(&scal_global(1,1),ngl*NumEdges,MPI_DOUBLE,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);
    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);



ElementToEdge.resize(4,NumElements);
ElemEdgeMasterSlave.resize(4,NumElements);

    MPI_Irecv(&ElementToEdge(1,1),4*NumElements,MPI_INT,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);
    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);
    MPI_Irecv(&ElemEdgeMasterSlave(1,1),4*NumElements,MPI_INT,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);
    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);
//    for (int i=1;i<=NumEdges;i++){
//        cout << "rank: " << MPI.rank << " Local Ele ID: " << i << "   Global ID: " << MyElementLocalToGlobal(i) << "\n" ;
//        for (int j=1;j<=10;j++){
//            cout << "rank: " << MPI.rank << " Local Edge ID " << i << "   EdgeInfo: "<<j <<" : " << MyEdgeInfo(j,i) << "\n" ;
//        }
//
//    }

}






void MeshPartitioning::SortMPIEdges(const MPI_setup MPI){


//    cout <<"I AM " << MPI.rank << " My Global Edges are: ";
//    for (int is = 1; is<=NumEdges; is++){
//      int globalEdgeID  = 	MyEdgeInfo(10,is);
//      cout << " " << globalEdgeID ;
//    }
//      cout <<"\n";

    MPIEdges.resize(2,NumProcessors);
    for (int i=1;i<=NumProcessors;i++){
        MPIEdges(1,i) = 0;
        MPIEdges(2,i) = i-1;

    }



  for (int is = 1; is<=NumEdges; is++){
//    cout << "Ele Local: "<< ie+1 << " Ele Global: "  << MeshSplit.ElementLocalToGlobal(ie+1,1) << "\n";
      int cpuL	        =	MyEdgeInfo(8,is);
      int cpuR	        =   MyEdgeInfo(9,is);
      int globalEdgeID  = 	MyEdgeInfo(10,is);
      // CHECK IF THIS IS A MPI EDGE!
      if ((cpuL != cpuR) && (cpuL!=-1) && (cpuR!=-1)){
        if (MPI.rank == cpuL){
            MPIEdges(1,cpuR+1) ++;
            MPIEdges(2,cpuR+1) = cpuR;

        }else{
            MPIEdges(1,cpuL+1) ++;
            MPIEdges(2,cpuL+1) = cpuL;

        }
      }else{
            MPIEdges(1,MPI.rank+1) ++;
            MPIEdges(2,MPI.rank+1) =MPI.rank;
      }
  }




  ProcIndex.resize(2,NumProcessors);
  ProcIndex(1,1) = 0;
  ProcIndex(2,NumProcessors) = NumEdges-1;
    for (int i=1;i<NumProcessors;i++){
        ProcIndex(2,i) = ProcIndex(1,i) + MPIEdges(1,i)-1;
        ProcIndex(1,i+1) = ProcIndex(2,i) + 1;

    }

//  for (int i=1;i<=NumProcessors;i++){
//    cout << "I am Rank "<<MPI.rank <<" and i have "<< MPIEdges(1,i) <<  " MPI Edges with Processor " <<MPIEdges(2,i) <<" storing edges from "<< ProcIndex(1,i) <<" to "<< ProcIndex(2,i) <<"\n";
//
//  }

    imatrix EdgeCounterPerProc(NumProcessors);
    for (int i=1;i<=NumProcessors;i++){
        EdgeCounterPerProc(i)= 0;
    }

    imatrix newEdgeIndex(NumEdges);
    imatrix oldEdgeIndex(NumEdges);
    imatrix EdgeInfoTMP(10,NumEdges);
    imatrix MyEdgesLocalToGlobalTMP(NumEdges);
    fmatrix nx_globalTMP(ngl*NumEdges);
    fmatrix ny_globalTMP(ngl*NumEdges);
    fmatrix scal_globalTMP(ngl*NumEdges);
    imatrix ElementToEdgeTMP(4,NumElements);
    imatrix ElemEdgeMasterSlaveTMP(4,NumElements);
    imatrix SwitchLeftRight(NumEdges);



    int sortCPU;
    for (int is = 1; is<=NumEdges; is++){
    //    cout << "Ele Local: "<< ie+1 << " Ele Global: "  << MeshSplit.ElementLocalToGlobal(ie+1,1) << "\n";
        int cpuL	        =	MyEdgeInfo(8,is);
        int cpuR	        =   MyEdgeInfo(9,is);
        int globalEdgeID  = 	MyEdgeInfo(10,is);
        SwitchLeftRight(is) = 0;



        if (cpuL==-1){sortCPU = cpuR+1;SwitchLeftRight(is)=1;}
        else if (cpuR==-1){sortCPU = cpuL+1;}
        else if (MPI.rank ==cpuL){sortCPU=cpuR+1;}
        else { sortCPU=cpuL+1;SwitchLeftRight(is)=1;}       // WE ARE RIGHT CPU! CHANGE THAT!!




        newEdgeIndex(is) =ProcIndex(1,sortCPU)+EdgeCounterPerProc(sortCPU)+1;
        oldEdgeIndex(newEdgeIndex(is))= is;
        EdgeCounterPerProc(sortCPU)++;
    }

//        cout << "RESTORING: old index: "<< is << " double check: " << oldEdgeIndex(newEdgeIndex(is)) << " and the new edge index is " <<newEdgeIndex(is)  <<"\n";
      //RESORT EDGE INFORMATION
//      for (int info = 1; info<=10;info++){
//            EdgeInfoTMP(info,newEdgeIndex(is)) = MyEdgeInfo(info,is);
//        }




int newOrderGlobal[NumEdges];
for (int iproc = 1; iproc <= NumProcessors;iproc++){
    int NumberEdges = ProcIndex(2,iproc)-ProcIndex(1,iproc)+1;


    for (int j=0;j<NumberEdges;j++){
        newOrderGlobal[ProcIndex(1,iproc)+j] = MyEdgeInfo(10,oldEdgeIndex(ProcIndex(1,iproc)+j+1));
    }
    sort(newOrderGlobal + ProcIndex(1,iproc),newOrderGlobal + ProcIndex(1,iproc)+NumberEdges);

}

//    cout <<"RANK: " << MPI.rank << " Edge Numbering  NEW : " ;
//    for (int j=0;j<NumEdges;j++){
//        cout << newOrderGlobal[j] << " " ;
//    }
//    cout <<"\n";
//
//        cout <<"RANK: " << MPI.rank << " Edge Numbering OLD  : " ;
//    for (int j=1;j<=NumEdges;j++){
//        cout << MyEdgeInfo(10,oldEdgeIndex(j)) << " " ;
//    }
//    cout <<"\n";



// RESORT newEdgeIndex array via the found global edge ordering on this processor
for (int is = 1; is<=NumEdges; is++){
    for (int j =0; j<NumEdges;j++){
        if (MyEdgeInfo(10,is) == newOrderGlobal[j]){
            newEdgeIndex(is) = j+1;
        }
    }
}


//    cout <<"RANK: " << MPI.rank << " Edge Numbering  NEW : " ;
//    for (int j=1;j<=NumEdges;j++){
//        cout << newEdgeIndex(j) << " " ;
//    }
//    cout <<"\n";






    for (int is = 1; is<=NumEdges; is++){
      EdgeInfoTMP(1,newEdgeIndex(is)) = MyEdgeInfo(1,is);
      EdgeInfoTMP(2,newEdgeIndex(is)) = MyEdgeInfo(2,is);
      if (SwitchLeftRight(is)==0){
          for (int info = 3; info<=10;info++){
                EdgeInfoTMP(info,newEdgeIndex(is)) = MyEdgeInfo(info,is);
            }

      }else{
    //            cout << "we switch left right!\n";
            EdgeInfoTMP(3,newEdgeIndex(is)) = MyEdgeInfo(4,is);
            EdgeInfoTMP(4,newEdgeIndex(is)) = MyEdgeInfo(3,is);
            EdgeInfoTMP(5,newEdgeIndex(is)) = MyEdgeInfo(6,is);
            EdgeInfoTMP(6,newEdgeIndex(is)) = MyEdgeInfo(5,is);
            EdgeInfoTMP(7,newEdgeIndex(is)) = MyEdgeInfo(7,is);
            EdgeInfoTMP(8,newEdgeIndex(is)) = MyEdgeInfo(9,is);
            EdgeInfoTMP(9,newEdgeIndex(is)) = MyEdgeInfo(8,is);
            EdgeInfoTMP(10,newEdgeIndex(is)) = MyEdgeInfo(10,is);
      }

      MyEdgesLocalToGlobalTMP(newEdgeIndex(is)) = EdgeInfoTMP(10,newEdgeIndex(is));


      for (int i = 1;i<=ngl;i++){
        if (SwitchLeftRight(is) ==0){
            nx_globalTMP((newEdgeIndex(is)-1)*ngl+i) = nx_global((is-1)*ngl+i);
            ny_globalTMP((newEdgeIndex(is)-1)*ngl+i) = ny_global((is-1)*ngl+i);
        }else{
            nx_globalTMP((newEdgeIndex(is)-1)*ngl+i) = -nx_global((is-1)*ngl+i);
            ny_globalTMP((newEdgeIndex(is)-1)*ngl+i) = -ny_global((is-1)*ngl+i);
        }

        scal_globalTMP((newEdgeIndex(is)-1)*ngl+i) = scal_global((is-1)*ngl+i);
      }

  }

  for (int ie=1;ie<=NumElements;ie++){
    for (int is=1;is<=4;is++){
        ElementToEdgeTMP(is,ie) = newEdgeIndex(ElementToEdge(is,ie));
        if (SwitchLeftRight(ElementToEdge(is,ie)) ==0){
            ElemEdgeMasterSlaveTMP(is,ie) = ElemEdgeMasterSlave(is,ie);
        }else{
            ElemEdgeMasterSlaveTMP(is,ie) = -ElemEdgeMasterSlave(is,ie);
        }
    }

  }
//    if (MPI.rank==0){
//    for (int is = 1; is<=NumEdges; is++){
////    cout << "Ele Local: "<< ie+1 << " Ele Global: "  << MeshSplit.ElementLocalToGlobal(ie+1,1) << "\n";
//      int cpuL	        =	MyEdgeInfo(8,is);
//      int cpuR	        =   MyEdgeInfo(9,is);
//      int globalEdgeID  = 	MyEdgeInfo(10,is);
//
//    cout << "global Edge id: "  << globalEdgeID << " cpuL: " << cpuL << " cpuR: " << cpuR << "\n";
//
//    }}
    MyEdgeInfo=EdgeInfoTMP;
    nx_global=nx_globalTMP;
    ny_global=ny_globalTMP;
    scal_global=scal_globalTMP;
    ElementToEdge = ElementToEdgeTMP;
    ElemEdgeMasterSlave= ElemEdgeMasterSlaveTMP;
    MyEdgesLocalToGlobal= MyEdgesLocalToGlobalTMP;



//   for (int ie=1;ie<=NumElements;ie++){
//        cout << " I am Rank: " << MPI.rank << " this is Element: " << ie ;
//    for (int is=1;is<=4;is++){
//        cout << " side " << ElementToEdgeTMP(is,ie) << " SwitchLeftRight is " << SwitchLeftRight(ElementToEdge(is,ie)) <<"   ";
//
//    }
//    cout <<"\n";
//  }


//        if (MPI.rank==0){
//        for (int is = 1; is<=NumEdges; is++){
////    cout << "Ele Local: "<< ie+1 << " Ele Global: "  << MeshSplit.ElementLocalToGlobal(ie+1,1) << "\n";
//      int cpuL	        =	MyEdgeInfo(8,is);
//      int cpuR	        =   MyEdgeInfo(9,is);
//      int globalEdgeID  = 	MyEdgeInfo(10,is);
//
//    cout << "global Edge id: "  << globalEdgeID << " cpuL: " << cpuL << " cpuR: " << cpuR << "\n";
//
//    }}


//    for (int is = 1; is<=NumEdges; is++){
//        cout << "My local Edge " << is << " is global edge "<< EdgeInfoTMP(10,is) << " and used to be local edge "<< oldEdgeIndex(is) <<"\n" ;
//    }




//    for (int i=1;i<=NumProcessors;i++){
//    cout << "I am Rank "<<MPI.rank <<" and i have "<< EdgeCounterPerProc(i) <<  " Edges with Processor " <<i <<" storing edges from "<< ProcIndex(1,i) <<" to "<< ProcIndex(2,i) <<"\n";
//
//  }


}














void MeshPartitioning::isInArray(const int SearchIndex, const int IntArray[],const int ArraySize, int *FoundLocation){

*FoundLocation = 0;
for (int i=1;i<=ArraySize;i++){
    if(SearchIndex == IntArray[i-1]){
        *FoundLocation=i;
    }
}
}




// Find some approximation to each element size
void MeshPartitioning :: ApproximateElementSizes(const int NumElements,const int ngl2, const dfloat x_phy[], const dfloat y_phy[], dfloat ElementSizes[]){

for (int ie=0;ie<NumElements;ie++){
//IDEA: find shortest edge and approximate element size as a square with that edge length
int id = ie*ngl2;
//Eck1
dfloat x1 = x_phy[id];
dfloat y1 = y_phy[id];
//Eck2
dfloat x2 = x_phy[id+ngl-1];
dfloat y2 = y_phy[id+ngl-1];
//Eck3
dfloat x3 = x_phy[id+(ngl-1)*ngl+ngl-1];
dfloat y3 = y_phy[id+(ngl-1)*ngl+ngl-1];
//Eck3
dfloat x4 = x_phy[id+(ngl-1)*ngl];
dfloat y4 = y_phy[id+(ngl-1)*ngl];


dfloat LengthEdge1 = sqrt(pow((x2-x1),2) + pow((y2-y1),2));
dfloat LengthEdge2 = sqrt(pow((x3-x2),2) + pow((y3-y2),2));
dfloat LengthEdge3 = sqrt(pow((x4-x3),2) + pow((y4-y3),2));
dfloat LengthEdge4 = sqrt(pow((x4-x1),2) + pow((y4-y1),2));

dfloat MinLength = min(min(LengthEdge1,LengthEdge2),min(LengthEdge3,LengthEdge4));

ElementSizes[ie] = MinLength;//*MinLength;

}




}

