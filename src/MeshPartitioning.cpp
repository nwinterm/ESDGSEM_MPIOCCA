#include "MeshPartitioning.h"


using namespace std;
MeshPartitioning::MeshPartitioning( const int procs, const int N)
{
    NumProcessors=procs;
    ngl = N;// N passed here is really N+1
    ngl2 = N*N;
    CommTags = (int*) calloc(NumProcessors*NumProcessors,sizeof(int)); //.resize(NumProcessors,NumProcessors);
    int counter=0;
    for (int i=0; i<NumProcessors; i++)
    {
        for (int j=0; j<NumProcessors; j++)
        {
            CommTags[i*NumProcessors + j ] = counter;
            counter++;
        }
    }

}

MeshPartitioning::~MeshPartitioning()
{
    //dtor
}

void MeshPartitioning::DivideBottom(const MPI_setup MPI, const dfloat b_global, dfloat * b)
{
    if (NumProcessors==1)
    {
        std::memcpy(&b, &b_global, sizeof b_global);
    }
    else
    {
        SplitRealValuesBetweenProcs( MPI,b_global, b);
    }


}


void MeshPartitioning::DivideMesh(const Mesh GlobalMesh,const MPI_setup MPI)
{


    if (NumProcessors==1)
    {



        global_NumElements=GlobalMesh.m_num_elements;
        global_NumEdges=GlobalMesh.m_num_edges;
        NumElements = global_NumElements;
        NumEdges=global_NumEdges;

        ElementsPerProc= (int*) calloc(1,sizeof(int));
        EdgesPerProc = (int*) calloc(1,sizeof(int));

        ElementsPerProc[0] = global_NumElements;
        EdgesPerProc[0] = global_NumEdges;

        ElementLocalToGlobal= (int*) calloc(ElementsPerProc[0]*1,sizeof(int));

        for (int i = 0; i<NumElements; i++)
        {
            ElementLocalToGlobal[i]=i+1;
        }


        cout <<"using " << NumProcessors <<" processors for a problem with "<< global_NumElements <<" Elements. and "<< global_NumEdges << "edges. \n";

        MyElementLocalToGlobal = (int*) calloc(NumElements,sizeof(int));
        for( int ie = 0; ie<NumElements; ie++)
        {
            MyElementLocalToGlobal[ie] = ie+1;
        }




        MyEdgesLocalToGlobal= (int*) calloc(NumEdges,sizeof(int));

        MyEdgeInfo = (int*) calloc(10*NumEdges,sizeof(int));

        for (int is=0; is<NumEdges; is++)
        {
            MyEdgesLocalToGlobal[is] = is; //(is,1);
        }

        for (int is=0; is<global_NumEdges; is++)
        {
            int glbEdgeInfoID=(is)*7;
            int id1 =  10*(is);
            MyEdgeInfo[id1+0]= GlobalMesh.EdgeInfo[glbEdgeInfoID+0]	;		//	!these are not even used after the inital mesh generation
            MyEdgeInfo[id1+1] = GlobalMesh.EdgeInfo[glbEdgeInfoID+1]	;		//	! ""
            MyEdgeInfo[id1+2] = GlobalMesh.EdgeInfo[glbEdgeInfoID+2] ;//!this is changed to be the local element on left processor
            MyEdgeInfo[id1+3] = GlobalMesh.EdgeInfo[glbEdgeInfoID+3] ;//!this is changed to be the local element on right processor
            MyEdgeInfo[id1+4] = GlobalMesh.EdgeInfo[glbEdgeInfoID+4];	//!local side left element <- fine as is
            MyEdgeInfo[id1+5] = GlobalMesh.EdgeInfo[glbEdgeInfoID+5]	;//!local side right element <- fine as is
            MyEdgeInfo[id1+6] = GlobalMesh.EdgeInfo[glbEdgeInfoID+6];

            MyEdgeInfo[id1+7]= 0;
            if (MyEdgeInfo[id1+3]==-1)
            {
                MyEdgeInfo[id1+8]= -1;
            }
            else
            {
                MyEdgeInfo[id1+8]= 0;
            }

            MyEdgeInfo[id1+9] = is 	;	//!global side id
        }




//        std::memcpy(&MyEdgeInfo, &GlobalMesh.EdgeInfo, sizeof GlobalMesh.EdgeInfo);





        x_global = (dfloat*) calloc(NumElements*ngl2,sizeof(dfloat));
        y_global = (dfloat*) calloc(NumElements*ngl2,sizeof(dfloat));
        xXi_global = (dfloat*) calloc(NumElements*ngl2,sizeof(dfloat));
        xEta_global = (dfloat*) calloc(NumElements*ngl2,sizeof(dfloat));
        yXi_global = (dfloat*) calloc(NumElements*ngl2,sizeof(dfloat));
        yEta_global = (dfloat*) calloc(NumElements*ngl2,sizeof(dfloat));
        J_global = (dfloat*) calloc(NumElements*ngl2,sizeof(dfloat));


        std::memcpy(&x_global, &GlobalMesh.x_global, sizeof GlobalMesh.x_global);
        std::memcpy(&y_global, &GlobalMesh.y_global, sizeof GlobalMesh.y_global);
        std::memcpy(&xXi_global, &GlobalMesh.xXi_global, sizeof GlobalMesh.xXi_global);
        std::memcpy(&xEta_global, &GlobalMesh.xEta_global, sizeof GlobalMesh.xEta_global);
        std::memcpy(&yXi_global, &GlobalMesh.yXi_global, sizeof GlobalMesh.yXi_global);
        std::memcpy(&yEta_global, &GlobalMesh.yEta_global, sizeof GlobalMesh.yEta_global);
        std::memcpy(&J_global, &GlobalMesh.J_global, sizeof GlobalMesh.J_global);



        nx_global = (dfloat*) calloc(ngl*NumEdges,sizeof(dfloat));
        ny_global = (dfloat*) calloc(ngl*NumEdges,sizeof(dfloat));
        scal_global = (dfloat*) calloc(ngl*NumEdges,sizeof(dfloat));

        std::memcpy(&nx_global, &GlobalMesh.NormalsX, sizeof GlobalMesh.NormalsX);
        std::memcpy(&ny_global, &GlobalMesh.NormalsY, sizeof GlobalMesh.NormalsY);
        std::memcpy(&scal_global, &GlobalMesh.Scal, sizeof GlobalMesh.Scal);




        MyElementToEdge = (int*) calloc(4*NumElements,sizeof(int));    //.resize(4,ElementsPerProc[0]);
        MyElemEdgeMasterSlave= (int*) calloc(4*NumElements,sizeof(int));    //.resize(4,ElementsPerProc[0]);

        std::memcpy(&MyElementToEdge, &GlobalMesh.ElementToEdge, sizeof GlobalMesh.ElementToEdge);
        std::memcpy(&MyElemEdgeMasterSlave, &GlobalMesh.ElemEdgeMasterSlave, sizeof GlobalMesh.ElemEdgeMasterSlave);




    }
    else
    {


        global_NumElements=GlobalMesh.m_num_elements;

        MPI_Bcast(&global_NumElements,1,MPI_INT,0,MPI_COMM_WORLD);
        global_NumEdges=GlobalMesh.m_num_edges;

        ElementsPerProc= (int*) calloc(NumProcessors,sizeof(int));
        EdgesPerProc = (int*) calloc(NumProcessors,sizeof(int));


        cout <<"using " << NumProcessors <<" processors for a problem with "<< global_NumElements <<" Elements.\n";
        for (int i=0; i<NumProcessors; i++)
        {
            ElementsPerProc[i] = global_NumElements / NumProcessors;
        }

        int RestElements = global_NumElements % NumProcessors;

        for (int i=0; i<RestElements; i++)
        {
            ElementsPerProc[i]=ElementsPerProc[i]+1;
        }


        MPI_Scatter (&ElementsPerProc[0],1,MPI_INT,&NumElements,1,MPI_INT,0,MPI_COMM_WORLD);

        ElementLocalToGlobal= (int*) calloc(ElementsPerProc[0]*NumProcessors,sizeof(int));
        //.resize(ElementsPerProc[0],NumProcessors);  //Link LOCAL Elements for EACH Processors to GLOBAL ELEMENT ID
        int * ElementGlobalToLocal = (int*) calloc(2*global_NumElements,sizeof(int));
        //.resize(2,global_NumElements); //Link each GLOBAL element to LOCALELEMENT + PROCESSOR

        int eleID_global=1;
        for (int iproc=1; iproc<=NumProcessors; iproc++)
        {
            for (int ie=1; ie<=ElementsPerProc[iproc-1]; ie++)
            {

                ElementLocalToGlobal[(iproc-1)*ElementsPerProc[0] + ie-1]= eleID_global; //(ie,iproc)
                ElementGlobalToLocal[eleID_global-1]=ie;//(1,eleID_global) = ie;
                ElementGlobalToLocal[global_NumElements + eleID_global-1] = iproc;
                eleID_global++;

            }
        }

        MyElementLocalToGlobal = (int*) calloc(NumElements,sizeof(int));
        for( int ie = 0; ie<ElementsPerProc[0]; ie++)
        {
            MyElementLocalToGlobal[ie] = ElementLocalToGlobal[ie];
        }

        for (int iproc=1; iproc<NumProcessors; iproc++)
        {
            MPI_Send (&ElementLocalToGlobal[iproc*ElementsPerProc[0] ],ElementsPerProc[iproc],MPI_INT,iproc,iproc,MPI_COMM_WORLD);
        }


        int * EdgeArray_proc = (int*) calloc(global_NumEdges,sizeof(int));

//int EdgeArray_proc[global_NumEdges];
        int LocatedAt;
        int maxEdgesLocal=0;
        EdgeLocalToGlobal= (int*) calloc(global_NumEdges*NumProcessors,sizeof(int));  //.resize(global_NumEdges,NumProcessors);
        int * EdgeGlobalToLocal= (int*) calloc(global_NumEdges*4,sizeof(int));  //.resize(4,global_NumEdges);

        for (int i=1; i<=global_NumEdges; i++)
        {
            EdgeArray_proc[i-1]=-1;
            for (int iproc=1; iproc<=NumProcessors; iproc++)
            {
                EdgeLocalToGlobal[(iproc-1)*global_NumEdges + i-1]=-1;   //(i,iproc)
            }
        }
        for (int iproc=0; iproc<NumProcessors; iproc++)
        {
            EdgesPerProc[iproc] = 0;
        }

        for (int iproc=1; iproc<=NumProcessors; iproc++)
        {

            // use an array to store all global edges that belong to this processor
            for (int i=1; i<=global_NumEdges; i++)
            {
                EdgeArray_proc[i-1] = -1; //EdgeLocalToGlobal(i,iproc);
            }

            for (int ie=1; ie<=ElementsPerProc[iproc-1]; ie++)
            {
                eleID_global = ElementLocalToGlobal[(iproc-1)*ElementsPerProc[0] + ie-1];

                for (int is=1; is<=4; is++)
                {
                    //get global side ID
                    int id = (eleID_global-1)*4 + is-1;
                    int sideID_global=GlobalMesh.ElementToEdge[id];





                    if (isInArray2(sideID_global,EdgeArray_proc,EdgesPerProc[iproc-1]) == 1)
                    {

                        if (GlobalMesh.EdgeInfo[(sideID_global-1)*7+2]==eleID_global-1)
                        {
//                        cout << " WE ACTUALLY GET IN WEIRD CASE WHERE RIGHT ELE WAS FOUND FIRST\n";
                            EdgeGlobalToLocal[(sideID_global-1)*4+0] = EdgeGlobalToLocal[(sideID_global-1)*4+1];
                            EdgeGlobalToLocal[(sideID_global-1)*4+2] = EdgeGlobalToLocal[(sideID_global-1)*4+3];
                        }
                        else
                        {
//                        if (sideID_global>3280){cout <<"STIMMT WAS NICHT!";}
                            EdgeGlobalToLocal[(sideID_global-1)*4+1] = EdgeGlobalToLocal[(sideID_global-1)*4];
                            EdgeGlobalToLocal[(sideID_global-1)*4+3] = EdgeGlobalToLocal[(sideID_global-1)*4+2];
                        }

                    }
                    else
                    {

                        //PRINT*, 'SIDE NOT FOUND: ', sideID_global
                        EdgesPerProc[iproc-1] = EdgesPerProc[iproc-1]+1;
//                    cout << " EDGE ADDED : " << sideID_global << "Element: " << eleID_global << " side " << is <<"\n" ;

                        EdgeLocalToGlobal[(iproc-1)*global_NumEdges + EdgesPerProc[iproc-1]-1] = sideID_global; //(EdgesPerProc[iproc-1],iproc)
                        EdgeArray_proc[EdgesPerProc[iproc-1]-1] = sideID_global; //EdgeLocalToGlobal(i,iproc);
                        //check if this element is LEFT or RIGHT element to the edge
                        if (GlobalMesh.EdgeInfo[(sideID_global-1)*7+2]==eleID_global-1)
                        {
                            EdgeGlobalToLocal[(sideID_global-1)*4] = EdgesPerProc[iproc-1];    //(1,sideID_global)
                            EdgeGlobalToLocal[(sideID_global-1)*4+2] = iproc;
                        }
                        else
                        {

                            EdgeGlobalToLocal[(sideID_global-1)*4+1] = EdgesPerProc[iproc-1];
                            EdgeGlobalToLocal[(sideID_global-1)*4+3] = iproc;
                        }
                    }

                }

            }
            maxEdgesLocal = max(maxEdgesLocal,EdgesPerProc[iproc-1])	;
        }


        int * globalEdgeInfo = (int*) calloc(10*maxEdgesLocal*NumProcessors,sizeof(int));

        int localEdgeCounter[NumProcessors];

        for (int iproc=0; iproc<NumProcessors; iproc++)
        {
            localEdgeCounter[iproc] = 0;
        }



        for (int is=1; is<=global_NumEdges; is++)
        {
            int edgeIDlocal_LEFT  = EdgeGlobalToLocal[(is-1)*4];		// 1=LOCAL SIDE ID PROC LEFT,  2=LOCAL SIDE ID PROC RIGHT,  3=LEFT PROCESSOR, 4=RIGHT PROCESSOR
            int edgeIDlocal_RIGHT = EdgeGlobalToLocal[(is-1)*4+1];
            int proc_LEFT	    = EdgeGlobalToLocal[(is-1)*4+2];
            int proc_RIGHT	    = EdgeGlobalToLocal[(is-1)*4+3];
            localEdgeCounter[proc_LEFT-1]++;
            int glbEdgeInfoID=(is-1)*7;

//CHECK IF THIS EDGE CROSSES PROCESSOR LINES
            if (proc_LEFT != proc_RIGHT)
            {
                if (proc_RIGHT == 0 )
                {
                    int id1 = 10*maxEdgesLocal*(proc_LEFT-1)+ 10*(edgeIDlocal_LEFT-1);

                    globalEdgeInfo[id1+0] = GlobalMesh.EdgeInfo[glbEdgeInfoID];//				!these are not even used after the inital mesh generation
                    globalEdgeInfo[id1+1] = GlobalMesh.EdgeInfo[glbEdgeInfoID+1];//			! ""
                    globalEdgeInfo[id1+2] = ElementGlobalToLocal[GlobalMesh.EdgeInfo[glbEdgeInfoID+2]] -1;//(1,GlobalMesh.EdgeInfo[glbEdgeInfoID+2]+1)-1	;//this is changed to be the local element on left processor
                    globalEdgeInfo[id1+3] = -1; //neighbour element is zero since its a boundary
                    globalEdgeInfo[id1+4] = GlobalMesh.EdgeInfo[glbEdgeInfoID+4]	;//local side left element <- fine as is
                    globalEdgeInfo[id1+5] = GlobalMesh.EdgeInfo[glbEdgeInfoID+5];	//local side right element <- fine as is
                    globalEdgeInfo[id1+6] = GlobalMesh.EdgeInfo[glbEdgeInfoID+6];
                    globalEdgeInfo[id1+7] = proc_LEFT -1;
                    globalEdgeInfo[id1+8] = proc_RIGHT -1;
                    globalEdgeInfo[id1+9] = is -1	;	//global side id
                }
                else
                {
                    localEdgeCounter[proc_RIGHT-1]++;

                    int id1 = 10*maxEdgesLocal*(proc_LEFT-1)+ 10*(edgeIDlocal_LEFT-1);
                    int id2 = 10*maxEdgesLocal*(proc_RIGHT-1)+ 10*(edgeIDlocal_RIGHT-1);

                    globalEdgeInfo[id1+0] = GlobalMesh.EdgeInfo[glbEdgeInfoID];				//these are not even used after the inital mesh generation
                    globalEdgeInfo[id1+1] = GlobalMesh.EdgeInfo[glbEdgeInfoID+1];		//! ""
                    globalEdgeInfo[id1+2] = ElementGlobalToLocal[GlobalMesh.EdgeInfo[glbEdgeInfoID+2]]-1;//(1,GlobalMesh.EdgeInfo[glbEdgeInfoID+2]+1)-1	;//!this is changed to be the local element on left processor
                    globalEdgeInfo[id1+3] = GlobalMesh.EdgeInfo[glbEdgeInfoID+3] 	;//!this is changed to be the GLOBAL ELEMENT ID of the right element!!!!
                    globalEdgeInfo[id1+4] = GlobalMesh.EdgeInfo[glbEdgeInfoID+4]	;//!local side left element <- fine as is
                    globalEdgeInfo[id1+5] = GlobalMesh.EdgeInfo[glbEdgeInfoID+5]	;//!local side right element <- fine as is
                    globalEdgeInfo[id1+6] = GlobalMesh.EdgeInfo[glbEdgeInfoID+6];
                    globalEdgeInfo[id1+7] = proc_LEFT -1;
                    globalEdgeInfo[id1+8] = proc_RIGHT -1;
                    globalEdgeInfo[id1+9] = is -1 	;	//!global side id

                    globalEdgeInfo[id2+0] = GlobalMesh.EdgeInfo[glbEdgeInfoID+0]	;	//		!these are not even used after the inital mesh generation
                    globalEdgeInfo[id2+1]  = GlobalMesh.EdgeInfo[glbEdgeInfoID+1];			//	! ""
                    globalEdgeInfo[id2+2]  = GlobalMesh.EdgeInfo[glbEdgeInfoID+2],	//!this is changed to be the GLOBAL ELEMENT ID on left processor
                                             globalEdgeInfo[id2+3]  = ElementGlobalToLocal[GlobalMesh.EdgeInfo[glbEdgeInfoID+3]]-1;//	!this is changed to be the local element on right processor
                    globalEdgeInfo[id2+4]  = GlobalMesh.EdgeInfo[glbEdgeInfoID+4];	//!local side left element <- fine as is
                    globalEdgeInfo[id2+5]  = GlobalMesh.EdgeInfo[glbEdgeInfoID+5];//!local side right element <- fine as is
                    // if orientation flip is on MPI edge, take special care!!
                    globalEdgeInfo[id2+6]  = GlobalMesh.EdgeInfo[glbEdgeInfoID+6];
                    //globalEdgeInfo[id2+6]  = GlobalMesh.EdgeInfo[glbEdgeInfoID+6];
                    globalEdgeInfo[id2+7]  = proc_LEFT -1;
                    globalEdgeInfo[id2+8]  = proc_RIGHT -1;
                    globalEdgeInfo[id2+9]  = is -1	;	//!global side id
                }
            }
            else
            {

                int id1 = 10*maxEdgesLocal*(proc_LEFT-1)+ 10*(edgeIDlocal_LEFT-1);
                globalEdgeInfo[id1+0]= GlobalMesh.EdgeInfo[glbEdgeInfoID+0]	;		//	!these are not even used after the inital mesh generation
                globalEdgeInfo[id1+1] = GlobalMesh.EdgeInfo[glbEdgeInfoID+1]	;		//	! ""
                globalEdgeInfo[id1+2] = ElementGlobalToLocal[GlobalMesh.EdgeInfo[glbEdgeInfoID+2] ]-1;//1,GlobalMesh.EdgeInfo[glbEdgeInfoID+2]+1)-1	; //!this is changed to be the local element on left processor
                globalEdgeInfo[id1+3] = ElementGlobalToLocal[GlobalMesh.EdgeInfo[glbEdgeInfoID+3] ]-1;//(1,GlobalMesh.EdgeInfo[glbEdgeInfoID+3]+1)-1 ;	//!this is changed to be the local element on right processor
                globalEdgeInfo[id1+4] = GlobalMesh.EdgeInfo[glbEdgeInfoID+4];	//!local side left element <- fine as is
                globalEdgeInfo[id1+5] = GlobalMesh.EdgeInfo[glbEdgeInfoID+5]	;//!local side right element <- fine as is
                globalEdgeInfo[id1+6] = GlobalMesh.EdgeInfo[glbEdgeInfoID+6];
                globalEdgeInfo[id1+7]= proc_LEFT -1;
                globalEdgeInfo[id1+8]= proc_LEFT -1;
                globalEdgeInfo[id1+9] = is -1	;	//!global side id

            }
        }



        MPI_Scatter (&EdgesPerProc[0],1,MPI_INT,&NumEdges,1,MPI_INT,0,MPI_COMM_WORLD);

        MyEdgesLocalToGlobal= (int*) calloc(NumEdges,sizeof(int));

        MyEdgeInfo = (int*) calloc(10*NumEdges,sizeof(int));

        for (int is=1; is<=EdgesPerProc[0]; is++)
        {

            MyEdgesLocalToGlobal[is-1] = EdgeLocalToGlobal[is-1]; //(is,1);

            int id =  10*(is-1);
            for (int j=0; j<10; j++)
            {
                MyEdgeInfo[id+j]= globalEdgeInfo[id+j];
            }

        }


        for (int iproc=1; iproc<NumProcessors; iproc++)
        {
            MPI_Send (&EdgeLocalToGlobal[iproc*global_NumEdges],EdgesPerProc[iproc],MPI_INT,iproc,iproc,MPI_COMM_WORLD);
        }
//    MPI_Waitall(NumProcessors,MPI.reqs, MPI.stats);

        for (int iproc=1; iproc<NumProcessors; iproc++)
        {

            int * globalEdgeInfo_TMP= (int*) calloc(10*EdgesPerProc[iproc],sizeof(int));
            for (int i=0; i<EdgesPerProc[iproc]; i++)
            {
                int id = 10*maxEdgesLocal*iproc+ 10*i;
                for (int j=0; j<10; j++)
                {
                    globalEdgeInfo_TMP[i*10+j]= globalEdgeInfo[id+j];
                }

            }
            MPI_Send (&globalEdgeInfo_TMP[0],10*EdgesPerProc[iproc],MPI_INT,iproc,iproc,MPI_COMM_WORLD);
//        MPI_Wait(&MPI.reqs[iproc], MPI.stats);
            free(globalEdgeInfo_TMP);
        }


        x_global = (dfloat*) calloc(NumElements*ngl2,sizeof(dfloat));
        y_global = (dfloat*) calloc(NumElements*ngl2,sizeof(dfloat));
        xXi_global = (dfloat*) calloc(NumElements*ngl2,sizeof(dfloat));
        xEta_global = (dfloat*) calloc(NumElements*ngl2,sizeof(dfloat));
        yXi_global = (dfloat*) calloc(NumElements*ngl2,sizeof(dfloat));
        yEta_global = (dfloat*) calloc(NumElements*ngl2,sizeof(dfloat));
        J_global = (dfloat*) calloc(NumElements*ngl2,sizeof(dfloat));



        SplitRealValuesBetweenProcs( MPI,GlobalMesh.x_global, x_global);

        SplitRealValuesBetweenProcs( MPI,GlobalMesh.y_global, y_global);

        SplitRealValuesBetweenProcs( MPI,GlobalMesh.yEta_global, yEta_global);

        SplitRealValuesBetweenProcs( MPI,GlobalMesh.xEta_global, xEta_global);

        SplitRealValuesBetweenProcs( MPI,GlobalMesh.yXi_global, yXi_global);

        SplitRealValuesBetweenProcs( MPI,GlobalMesh.xXi_global, xXi_global);

        SplitRealValuesBetweenProcs( MPI,GlobalMesh.J_global, J_global);



        nx_global = (dfloat*) calloc(ngl*NumEdges,sizeof(dfloat));
        ny_global = (dfloat*) calloc(ngl*NumEdges,sizeof(dfloat));
        scal_global = (dfloat*) calloc(ngl*NumEdges,sizeof(dfloat));

        SplitEdgeRealValuesBetweenProcs( MPI,GlobalMesh.NormalsX, nx_global);

        SplitEdgeRealValuesBetweenProcs( MPI,GlobalMesh.NormalsY, ny_global);

        SplitEdgeRealValuesBetweenProcs( MPI,GlobalMesh.Scal, scal_global);


        int * SplitElementToEdge = (int*) calloc(4*ElementsPerProc[0]*NumProcessors,sizeof(int));
        int * SplitElemEdgeMasterSlave = (int*) calloc(4*ElementsPerProc[0]*NumProcessors,sizeof(int));

        for (int ie=1; ie<=global_NumElements; ie++)
        {

            int eleID_local=ElementGlobalToLocal[(ie-1)];//(1,ie);
            int iproc= ElementGlobalToLocal[global_NumElements+(ie-1)];//(2,ie);

            for (int is=1; is<=4; is++)
            {
                int edgeID = GlobalMesh.ElementToEdge[(ie-1)*4+is-1];
                int id = 4*(iproc-1)*ElementsPerProc[0] + 4*(eleID_local-1) + (is-1);
                if (GlobalMesh.EdgeInfo[(edgeID-1)*7+2] == ie-1)
                {
                    SplitElementToEdge[id] = EdgeGlobalToLocal[(edgeID-1)*4];
//            SplitElementToEdge[is-1][eleID_local-1][iproc-1]  = EdgeGlobalToLocal(1,edgeID);
                }
                else
                {
                    SplitElementToEdge[id] = EdgeGlobalToLocal[(edgeID-1)*4+1];
//            SplitElementToEdge[is-1][eleID_local-1][iproc-1]  = EdgeGlobalToLocal(2,edgeID);
                }
                SplitElemEdgeMasterSlave[id] = GlobalMesh.ElemEdgeMasterSlave[(ie-1)*4+is-1];

            }
        }


        MyElementToEdge = (int*) calloc(4*ElementsPerProc[0],sizeof(int));    //.resize(4,ElementsPerProc[0]);
        MyElemEdgeMasterSlave= (int*) calloc(4*ElementsPerProc[0],sizeof(int));    //.resize(4,ElementsPerProc[0]);




        for (int ie=0; ie<ElementsPerProc[0]; ie++)
        {
            for (int is=0; is<4; is++)
            {
                int id =  4*ie + is;
                MyElementToEdge[id] = SplitElementToEdge[id];
                MyElemEdgeMasterSlave[id] = SplitElemEdgeMasterSlave[id];
            }
        }



        for (int iproc=2; iproc<=NumProcessors; iproc++)
        {
            int * EleToEdge_TMP = (int*) calloc(4*ElementsPerProc[iproc-1],sizeof(int));    // (4,ElementsPerProc[iproc-1]);
            int * EleToEdgeMS_TMP= (int*) calloc(4*ElementsPerProc[iproc-1],sizeof(int));    //(4,ElementsPerProc[iproc-1]);
            for (int ie=0; ie<ElementsPerProc[iproc-1]; ie++)
            {
                for (int is=0; is<4; is++)
                {
                    int id1 =  4*ie + is;
                    int id = 4*(iproc-1)*ElementsPerProc[0] + 4*ie + is;
                    EleToEdge_TMP[id1] = SplitElementToEdge[id];
//            EleToEdge_TMP(is,ie) = SplitElementToEdge[is-1][ie-1][iproc-1];
                    EleToEdgeMS_TMP[id1] = SplitElemEdgeMasterSlave[id];
                }
            }

            MPI_Send (&EleToEdge_TMP[0],4*ElementsPerProc[iproc-1],MPI_INT,iproc-1,iproc-1,MPI_COMM_WORLD);
            free(EleToEdge_TMP);

            MPI_Send (&EleToEdgeMS_TMP[0],4*ElementsPerProc[iproc-1],MPI_INT,iproc-1,iproc-1,MPI_COMM_WORLD);
            free(EleToEdgeMS_TMP);
        }

        free(EdgeArray_proc);
        free(globalEdgeInfo);
        free(SplitElementToEdge);
        free(SplitElemEdgeMasterSlave);

        free(ElementGlobalToLocal);
        free(EdgeGlobalToLocal);
    }

}

void MeshPartitioning::SplitRealValuesBetweenProcs(const MPI_setup MPI,const dfloat * Values, dfloat *HostValues)
{


    for (int ie=0; ie<ElementsPerProc[0]; ie++)
    {
        for (int i=0; i<ngl; i++)
        {
            for (int j=0; j<ngl; j++)
            {
                int id = (ElementLocalToGlobal[ie]-1)*ngl2   +j*ngl+i;
                int idLoc = ie*ngl2   +j*ngl+i;
                HostValues[idLoc] = Values[id];
            }
        }
    }





    for (int iproc=2; iproc<=NumProcessors; iproc++)
    {

        dfloat * tmpValues =    (dfloat*) calloc(ngl2*ElementsPerProc[iproc-1],sizeof(dfloat));
        for (int ie=0; ie<ElementsPerProc[iproc-1]; ie++)
        {
            for (int i=0; i<ngl; i++)
            {
                for (int j=0; j<ngl; j++)
                {
                    int id = (ElementLocalToGlobal[(iproc-1)*ElementsPerProc[0] + ie]-1)*ngl2   +j*ngl+i;
                    int idLoc = ie*ngl2   +j*ngl+i;
                    tmpValues[idLoc] = Values[id];
                }
            }
        }

        MPI_Send (&tmpValues[0],ngl2*ElementsPerProc[iproc-1],MPI_DFLOAT,iproc-1,iproc-1,MPI_COMM_WORLD);
//        MPI_Wait(&MPI.reqs[iproc-1], MPI.stats);
        free(tmpValues);

    }
}

void MeshPartitioning::SplitEdgeRealValuesBetweenProcs(const MPI_setup MPI,const dfloat * Values, dfloat * HostValues)
{


    for (int ie=0; ie<EdgesPerProc[0]; ie++)
    {
        for (int i=0; i<ngl; i++)
        {
            int id = (EdgeLocalToGlobal[ie]-1)*ngl   +i;
            int idLoc = ie*ngl  +i;
            HostValues[idLoc] = Values[id];
        }
    }





    for (int iproc=2; iproc<=NumProcessors; iproc++)
    {
        dfloat * tmpValues =    (dfloat*) calloc(ngl*EdgesPerProc[iproc-1],sizeof(dfloat));
        for (int ie=0; ie<EdgesPerProc[iproc-1]; ie++)
        {
            for (int i=0; i<ngl; i++)
            {
                int id = (EdgeLocalToGlobal[(iproc-1)*global_NumEdges +ie]-1)*ngl   +i;   // (ie+1,iproc)
                int idLoc = ie*ngl  +i;
                tmpValues[idLoc] = Values[id];
            }
        }

        MPI_Send (&tmpValues[0],ngl*EdgesPerProc[iproc-1],MPI_DFLOAT,iproc-1,iproc-1,MPI_COMM_WORLD);
//        MPI_Wait(&MPI.reqs[iproc-1], MPI.stats);
        free(tmpValues);

    }


}









































void MeshPartitioning::ReceiveBottom(const MPI_setup MPI, dfloat *b)
{

MPI_Recv(&b[0],ngl2*NumElements,MPI_DFLOAT,0,MPI.rank,MPI_COMM_WORLD, MPI.stats);

}




void MeshPartitioning::ReceiveMesh(const MPI_setup MPI)
{
    MPI_Bcast(&global_NumElements,1,MPI_INT,0,MPI_COMM_WORLD);

//    MPI_Recv(&global_NumElements,1,MPI_INT,0,MPI.rank,MPI_COMM_WORLD, MPI.stats);
//    MPI_Wait(&reqs[MPI.rank], stats);
//    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);


//    MPI_Recv(&NumElements,1,MPI_INT,0,MPI.rank,MPI_COMM_WORLD, MPI.stats);
    MPI_Scatter (&ElementsPerProc[0],1,MPI_INT,&NumElements,1,MPI_INT,0,MPI_COMM_WORLD);

//cout << "MPI rank: " << MPI.rank << " my elements: " << NumElements << "\n";
//    MPI_Wait(&reqs[MPI.rank], stats);
//    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);
//    MPI_Waitall(NumProcessors,MPI.reqs, MPI.stats);

    MyElementLocalToGlobal = (int*) calloc(NumElements,sizeof(int));


    MPI_Recv(&MyElementLocalToGlobal[0],NumElements,MPI_INT,0,MPI.rank,MPI_COMM_WORLD, MPI.stats);



    MPI_Scatter (&EdgesPerProc[0],1,MPI_INT,&NumEdges,1,MPI_INT,0,MPI_COMM_WORLD);


    MyEdgeInfo = (int*) calloc(10*NumEdges,sizeof(int));
    MyEdgesLocalToGlobal = (int*) calloc(NumEdges,sizeof(int));
    MPI_Recv(&MyEdgesLocalToGlobal[0],NumEdges,MPI_INT,0,MPI.rank,MPI_COMM_WORLD, MPI.stats);

    MPI_Recv(&MyEdgeInfo[0],10*NumEdges,MPI_INT,0,MPI.rank,MPI_COMM_WORLD, MPI.stats);




    x_global = (dfloat*) calloc(NumElements*ngl2,sizeof(dfloat));
    y_global = (dfloat*) calloc(NumElements*ngl2,sizeof(dfloat));
    xXi_global = (dfloat*) calloc(NumElements*ngl2,sizeof(dfloat));
    xEta_global = (dfloat*) calloc(NumElements*ngl2,sizeof(dfloat));
    yXi_global = (dfloat*) calloc(NumElements*ngl2,sizeof(dfloat));
    yEta_global = (dfloat*) calloc(NumElements*ngl2,sizeof(dfloat));
    J_global = (dfloat*) calloc(NumElements*ngl2,sizeof(dfloat));



    MPI_Recv(&x_global[0],ngl2*NumElements,MPI_DFLOAT,0,MPI.rank,MPI_COMM_WORLD, MPI.stats);
//    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);


    MPI_Recv(&y_global[0],ngl2*NumElements,MPI_DFLOAT,0,MPI.rank,MPI_COMM_WORLD, MPI.stats);
//    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);


    MPI_Recv(&yEta_global[0],ngl2*NumElements,MPI_DFLOAT,0,MPI.rank,MPI_COMM_WORLD, MPI.stats);
//    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);


    MPI_Recv(&xEta_global[0],ngl2*NumElements,MPI_DFLOAT,0,MPI.rank,MPI_COMM_WORLD, MPI.stats);
//    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);


    MPI_Recv(&yXi_global[0],ngl2*NumElements,MPI_DFLOAT,0,MPI.rank,MPI_COMM_WORLD, MPI.stats);
//    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);


    MPI_Recv(&xXi_global[0],ngl2*NumElements,MPI_DFLOAT,0,MPI.rank,MPI_COMM_WORLD, MPI.stats);
//    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);


    MPI_Recv(&J_global[0],ngl2*NumElements,MPI_DFLOAT,0,MPI.rank,MPI_COMM_WORLD, MPI.stats);
//    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);


    nx_global = (dfloat*) calloc(ngl*NumEdges,sizeof(dfloat));
    ny_global = (dfloat*) calloc(ngl*NumEdges,sizeof(dfloat));
    scal_global = (dfloat*) calloc(ngl*NumEdges,sizeof(dfloat));




    MPI_Recv(&nx_global[0],ngl*NumEdges,MPI_DFLOAT,0,MPI.rank,MPI_COMM_WORLD, MPI.stats);
//    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);

    MPI_Recv(&ny_global[0],ngl*NumEdges,MPI_DFLOAT,0,MPI.rank,MPI_COMM_WORLD, MPI.stats);
//    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);

    MPI_Recv(&scal_global[0],ngl*NumEdges,MPI_DFLOAT,0,MPI.rank,MPI_COMM_WORLD, MPI.stats);
//    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);



//ElementToEdge.resize(4,NumElements);
//ElemEdgeMasterSlave.resize(4,NumElements);

    MyElementToEdge = (int*) calloc(4*NumElements,sizeof(int));    //.resize(4,ElementsPerProc[0]);
    MyElemEdgeMasterSlave= (int*) calloc(4*NumElements,sizeof(int));    //.resize(4,ElementsPerProc[0]);

//    MPI_Irecv(&MyElementToEdge[0],4*NumElements,MPI_INT,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);
//    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);
//    MPI_Irecv(&MyElemEdgeMasterSlave[0],4*NumElements,MPI_INT,0,MPI.rank,MPI_COMM_WORLD,&MPI.reqs[MPI.rank]);
//    MPI_Wait(&MPI.reqs[MPI.rank], MPI.stats);

    MPI_Recv(&MyElementToEdge[0],4*NumElements,MPI_INT,0,MPI.rank,MPI_COMM_WORLD, MPI.stats);
    MPI_Recv(&MyElemEdgeMasterSlave[0],4*NumElements,MPI_INT,0,MPI.rank,MPI_COMM_WORLD, MPI.stats);
//    for (int i=0;i<NumElements;i++){
//        cout << "rank: " << MPI.rank << " Local Ele ID: " << i << "   edge 1: " << MyElementToEdge[i*4] << "   edge 2: " << MyElementToEdge[i*4+1] << "   edge 3: " << MyElementToEdge[i*4+2] << "   edge 2: " << MyElementToEdge[i*4+3] << "\n" ;
//
//
//    }

}






void MeshPartitioning::SortMPIEdges(const MPI_setup MPI)
{


    if (NumProcessors>1)
    {

//    MPIEdges.resize(2,NumProcessors);
        int * MPIEdges = (int*) calloc(2*NumProcessors,sizeof(int));
        // MPIEdges(1,:) = MPIEdges[0:NumProcessors-1]
        // MPIEdges(2,:) = MPIEdges[NumProcessors:2*NumProcessors-1]
        for (int i=1; i<=NumProcessors; i++)
        {
            MPIEdges[i-1] = 0;
            MPIEdges[NumProcessors+i-1] = i-1;

        }



        for (int is = 1; is<=NumEdges; is++)
        {
//    cout << "Ele Local: "<< ie+1 << " Ele Global: "  << MeshSplit.ElementLocalToGlobal(ie+1,1) << "\n";+
            int id = (is-1)*10;
            int cpuL	        =	MyEdgeInfo[id+7];
            int cpuR	        =   MyEdgeInfo[id+8];
            int globalEdgeID  = 	MyEdgeInfo[id+9];
            // CHECK IF THIS IS A MPI EDGE!
            if ((cpuL != cpuR) && (cpuL!=-1) && (cpuR!=-1))
            {
                if (MPI.rank == cpuL)
                {
                    MPIEdges[cpuR] ++;
                    MPIEdges[NumProcessors+cpuR] = cpuR;

                }
                else
                {
                    MPIEdges[cpuL] ++;
                    MPIEdges[NumProcessors+cpuL] = cpuL;

                }
            }
            else
            {
                MPIEdges[MPI.rank] ++;
                MPIEdges[NumProcessors+MPI.rank] =MPI.rank;
            }
        }


//  ProcIndex.resize(2,NumProcessors);
//  ProcIndex(1,1) = 0;
//  ProcIndex(2,NumProcessors) = NumEdges-1;
//    for (int i=1;i<NumProcessors;i++){
//        ProcIndex(2,i) = ProcIndex(1,i) + MPIEdges(1,i)-1;
//        ProcIndex(1,i+1) = ProcIndex(2,i) + 1;
//
//    }
//


        //ProcIndex stores the Edge IDs (first and last) for each processors
        ProcIndex = (int*) calloc(2*NumProcessors,sizeof(int));//.resize(2,NumProcessors);
        ProcIndex[0] = 0;
        ProcIndex[2*NumProcessors-1] = NumEdges-1;

        // 0 to NumProcessors-1  for start indices
        // NumProcessors to 2*NumProcessors-1 for end indices
        for (int i=1; i<NumProcessors; i++)
        {
            ProcIndex[NumProcessors + i-1] = ProcIndex[i-1] + MPIEdges[i-1]-1;
            ProcIndex[i] = ProcIndex[NumProcessors + i-1] + 1;

        }

//  for (int i=1;i<=NumProcessors;i++){
//    cout << "I am Rank "<<MPI.rank <<" and i have "<< MPIEdges(1,i) <<  " MPI Edges with Processor " <<MPIEdges(2,i) <<" storing edges from "<< ProcIndex[i-1] <<" to "<< ProcIndex[NumProcessors+i-1] <<"\n";
//
//  }


        int * EdgeInfoTMP= (int*) calloc(10*NumEdges,sizeof(int));

        int *  EdgeCounterPerProc = (int*) calloc(NumProcessors,sizeof(int));

        int *  newEdgeIndex = (int*) calloc(NumEdges,sizeof(int));
        int *  oldEdgeIndex = (int*) calloc(NumEdges,sizeof(int));
        int *  MyEdgesLocalToGlobalTMP = (int*) calloc(NumEdges,sizeof(int));
        int *  SwitchLeftRight = (int*) calloc(NumEdges,sizeof(int));


        dfloat * nx_globalTMP = (dfloat*) calloc(ngl*NumEdges,sizeof(dfloat));
        dfloat * ny_globalTMP = (dfloat*) calloc(ngl*NumEdges,sizeof(dfloat));
        dfloat * scal_globalTMP = (dfloat*) calloc(ngl*NumEdges,sizeof(dfloat));









        int sortCPU;
        for (int is = 1; is<=NumEdges; is++)
        {
            int id = (is-1)*10;
            //    cout << "Ele Local: "<< ie+1 << " Ele Global: "  << MeshSplit.ElementLocalToGlobal(ie+1,1) << "\n";
            int cpuL	        =	MyEdgeInfo[id+7];
            int cpuR	        =   MyEdgeInfo[id+8];
            int globalEdgeID  = 	MyEdgeInfo[id+9];
            SwitchLeftRight[is-1] = 0;



            if (cpuL==-1)
            {
                sortCPU = cpuR+1;
                SwitchLeftRight[is-1]=1;
            }
            else if (cpuR==-1)
            {
                sortCPU = cpuL+1;
            }
            else if (MPI.rank ==cpuL)
            {
                sortCPU=cpuR+1;
            }
            else
            {
                sortCPU=cpuL+1;    // WE ARE RIGHT CPU! CHANGE THAT!!
                SwitchLeftRight[is-1]=1;
            }




            newEdgeIndex[is-1] =ProcIndex[sortCPU-1]+EdgeCounterPerProc[sortCPU-1]+1;
            oldEdgeIndex[newEdgeIndex[is-1]-1]= is;
            EdgeCounterPerProc[sortCPU-1]++;
        }

        free(EdgeCounterPerProc);



        int newOrderGlobal[NumEdges];
        for (int iproc = 1; iproc <= NumProcessors; iproc++)
        {
            int NumberEdges = ProcIndex[NumProcessors+iproc-1]-ProcIndex[iproc-1]+1;


            for (int j=0; j<NumberEdges; j++)
            {
                int id = (oldEdgeIndex[ProcIndex[iproc-1]+j]-1)*10;
                newOrderGlobal[ProcIndex[iproc-1]+j] = MyEdgeInfo[id+9];
            }
            sort(newOrderGlobal + ProcIndex[iproc-1],newOrderGlobal + ProcIndex[iproc-1]+NumberEdges);

        }



// RESORT newEdgeIndex array via the found global edge ordering on this processor
        for (int is = 1; is<=NumEdges; is++)
        {
            for (int j =0; j<NumEdges; j++)
            {
                if (MyEdgeInfo[(is-1)*10+9] == newOrderGlobal[j])
                {
                    newEdgeIndex[is-1] = j+1;
                }
            }
        }






        for (int is = 1; is<=NumEdges; is++)
        {
            int id = (is-1)*10;
            int idTMP = (newEdgeIndex[is-1]-1)*10;
            EdgeInfoTMP[idTMP] = MyEdgeInfo[id];             //(1,newEdgeIndex[is-1])
            EdgeInfoTMP[idTMP+1] = MyEdgeInfo[id+1];
            if (SwitchLeftRight[is-1]==0)
            {
                for (int info = 2; info<10; info++)
                {
                    EdgeInfoTMP[idTMP+info]= MyEdgeInfo[id+info];
                }

            }
            else
            {
                //            cout << "we switch left right!\n";

                EdgeInfoTMP[idTMP+2] = MyEdgeInfo[id+3];
                EdgeInfoTMP[idTMP+3] = MyEdgeInfo[id+2];
                EdgeInfoTMP[idTMP+4] = MyEdgeInfo[id+5];
                EdgeInfoTMP[idTMP+5] = MyEdgeInfo[id+4];
                EdgeInfoTMP[idTMP+6] = MyEdgeInfo[id+6];
                EdgeInfoTMP[idTMP+7] = MyEdgeInfo[id+8];
                EdgeInfoTMP[idTMP+8] = MyEdgeInfo[id+7];
                EdgeInfoTMP[idTMP+9] = MyEdgeInfo[id+9];
            }

            MyEdgesLocalToGlobalTMP[newEdgeIndex[is-1]-1]= EdgeInfoTMP[idTMP+9];


            for (int i = 0; i<ngl; i++)
            {
                int oldID = (is-1)*ngl+i;
                int newID = (newEdgeIndex[is-1]-1)*ngl+i;
                if (SwitchLeftRight[is-1] ==0)
                {
                    nx_globalTMP[newID] = nx_global[oldID];
                    ny_globalTMP[newID] = ny_global[oldID];
                }
                else
                {
                    nx_globalTMP[newID] = -nx_global[oldID];
                    ny_globalTMP[newID] = -ny_global[oldID];
                }

                scal_globalTMP[newID] = scal_global[oldID];
            }

        }



        int * MyElementToEdgeTMP        = (int*) calloc(4*NumElements,sizeof(int));  //(4,NumElements);
        int * MyElemEdgeMasterSlaveTMP  = (int*) calloc(4*NumElements,sizeof(int));
        for (int ie=0; ie<NumElements; ie++)
        {
            for (int is=0; is<4; is++)
            {
                int id = ie*4 + is;
                MyElementToEdgeTMP[id] = newEdgeIndex[MyElementToEdge[id]-1];
                if (SwitchLeftRight[MyElementToEdge[id]-1] ==0)
                {
                    MyElemEdgeMasterSlaveTMP[id] = MyElemEdgeMasterSlave[id];
                }
                else
                {
                    MyElemEdgeMasterSlaveTMP[id] = -MyElemEdgeMasterSlave[id];
                }
            }
        }

        for (int is = 0; is<NumEdges; is++)
        {
            for (int info = 0; info<10; info++)
            {
                int id = is*10+info;
                MyEdgeInfo[id]=EdgeInfoTMP[id];
            }
            // TESTING BUGFIX FOR REVERSED ORIENTATION ON MPI EDGES!!!!!!!!!!!!!
//		int id = is*10+6;
//		if (MyEdgeInfo[id+2]!=MyEdgeInfo[id+1]){	// check if this is a MPI edge
//			if (MyEdgeInfo[is*10+6]==0 ){		// check if orientation is flipped
//                		if (SwitchLeftRight[is] ==1){	// check if left and right state where changed for this rank
//               				MyEdgeInfo[id]=-1;	// if this edge was Left/Right changed here and is an orientation flipped edge, note that now the left side is to be flipped!
//				}else{
//					MyEdgeInfo[id]=1;	// left/right state was not changed, but this is orientation flipped edge and MPI edge -> edge values will be flipped on OTHER rank and come to this rank correctly
//				}
//			}
//                }
        }

        for (int is = 0; is<NumEdges; is++)
        {
            for (int i = 0; i<ngl; i++)
            {
                int id= is*ngl+i;
                int id_new = id;

                if (SwitchLeftRight[is] ==1) 	// check if left and right state where changed for this rank
                {
                    if (MyEdgeInfo[is*10+6]==0 )
                    {
                        id_new = is*ngl+ngl-1-i;
                    }
                }


                nx_global[id_new]=nx_globalTMP[id];
                ny_global[id_new]=ny_globalTMP[id];
                scal_global[id_new]=scal_globalTMP[id];

            }
        }

        for (int ie=0; ie<NumElements; ie++)
        {
            for (int is=0; is<4; is++)
            {
                int id = ie*4 + is;
                MyElementToEdge[id]       = MyElementToEdgeTMP[id];
                MyElemEdgeMasterSlave[id] = MyElemEdgeMasterSlaveTMP[id];
            }
        }

        for (int is = 0; is<NumEdges; is++)
        {
            MyEdgesLocalToGlobal[is]= MyEdgesLocalToGlobalTMP[is];
        }


        free(nx_globalTMP);
        free(ny_globalTMP);
        free(scal_globalTMP);
        free(MPIEdges);
        free(MyElementToEdgeTMP);
        free(MyElemEdgeMasterSlaveTMP);
        free(newEdgeIndex);
        free(oldEdgeIndex);
        free(MyEdgesLocalToGlobalTMP);
        free(SwitchLeftRight);
        free(EdgeInfoTMP);
    }
}














void MeshPartitioning::isInArray(const int SearchIndex, const int IntArray[],const int ArraySize, int *FoundLocation)
{

    *FoundLocation = 0;
    for (int i=1; i<=ArraySize; i++)
    {
        if(SearchIndex == IntArray[i-1])
        {
            *FoundLocation=i;
        }
    }
}

int MeshPartitioning::isInArray2(const int SearchIndex, const int IntArray[],const int ArraySize)
{

    for (int i=1; i<=ArraySize; i++)
    {
        if(SearchIndex == IntArray[i-1])
        {
            return 1;
        }
    }
    return 0;
}

