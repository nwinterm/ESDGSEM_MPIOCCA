#include "Mesh.h"


Mesh::Mesh(const fmatrix fm_x_GL,const int int_ngl)
{
    x_GL=fm_x_GL;
    ngl=int_ngl;

}
Mesh::~Mesh()
{
    //dtor
}


void Mesh::InitMesh(const string meshFile, const bool Cartesian, const int Testcase)
{
NelemX=0;
NelemY=0;
if (Cartesian){
    dfloat xR,xL,yR,yL;
    bool PeriodicBD_X,PeriodicBD_Y;
    int fixedDomain,fixedDisc;


    InitDomain(Testcase, &fixedDomain,&fixedDisc,&NelemX, &NelemY,&PeriodicBD_X, &PeriodicBD_Y, &xL, &xR, &yL, &yR);
    ReadCartesianData(fixedDomain,fixedDisc,&xL,&xR,&yL,&yR,&NelemX,&NelemY,&PeriodicBD_X,&PeriodicBD_Y);
    GenerateMesh(xL,xR,yL,yR,PeriodicBD_X,PeriodicBD_Y);
    }
    else{
    ReadMesh(meshFile);
    }

    ElemEdgeMasterSlave.resize(4,m_num_elements);
    ElemEdgeOrientation.resize(4,m_num_elements);

    for(int ie=1;ie<=m_num_elements;++ie){
            for (int is=1;is<=4;is++){
                int ifa  = ElementToEdge(is,ie);

                if (EdgeInfo(3,ifa)==ie-1){
                    //this is the left element to this edge!
                    ElemEdgeMasterSlave(is,ie)=+1;
                    ElemEdgeOrientation(is,ie)=1;   //order is never reversed for left element! as it is master

                }else{
                    ElemEdgeMasterSlave(is,ie)=-1;
                    ElemEdgeOrientation(is,ie)=EdgeInfo(7,ifa);
                }

            }
      }





      // SPECIAL INTERIOR BOUNDARIES FOR EXAMPLE CURVED DAM BREAK




    if(Testcase ==43){   //PARTIAL CURVED DAM BREAK
        int ExtraEdges=36;
        int ghostCounter = 0;
        imatrix EdgeInfo_TMP(7,m_num_edges);

        for(int is=1;is<=m_num_edges;++is){
                for (int i=1;i<=7;i++){
                    EdgeInfo_TMP(i,is) = EdgeInfo(i,is);
                }
        }

        EdgeInfo.resize(7,m_num_edges+ExtraEdges);

        for(int is=1;is<=m_num_edges;++is){
                for (int i=1;i<=7;i++){
                      EdgeInfo(i,is) = EdgeInfo_TMP(i,is);
                }
        }
//
//
        for(int is=1;is<=m_num_edges;++is){

            if( (EdgeInfo(3,is) % 20 == 2) && (EdgeInfo(3,is) % 40 != 2) && (EdgeInfo(4,is) % 20 == 3)&& (EdgeInfo(4,is) % 40 != 3)){
                if ( (is!=1544) && (is!=1625) && (is != 1706) && (is!=1787) ){  //&& (is !=1868)

                ghostCounter = ghostCounter+1;
                EdgeInfo(1,m_num_edges+ghostCounter)= EdgeInfo(1,is);
                EdgeInfo(2,m_num_edges+ghostCounter)= EdgeInfo(2,is);
                EdgeInfo(3,m_num_edges+ghostCounter)= EdgeInfo(4,is);
                EdgeInfo(4,m_num_edges+ghostCounter)= -1;
                EdgeInfo(5,m_num_edges+ghostCounter)= EdgeInfo(6,is);
                EdgeInfo(6,m_num_edges+ghostCounter)= -1;
                EdgeInfo(7,m_num_edges+ghostCounter)= EdgeInfo(7,is);

                EdgeInfo(4,is) = -1 ;
                EdgeInfo(6,is)=-1;
                }
            }
        }
//        cout << "ghostcounter is " << ghostCounter <<"\n";
        ghostCounter = 1;
        for(int ie=24;ie<=1584;ie=ie+40){
            if ( (ie!=744) && (ie!=784) && (ie != 824) && (ie!=864) ){  //&& (ie !=903)
            ElementToEdge(4,ie)   = m_num_edges + ghostCounter;
            ElemEdgeMasterSlave(4,ie)=+1;
            ElemEdgeOrientation(4,ie)=1;   //order is never reversed for left element! as it is master

            ghostCounter = ghostCounter+1;
            }
        }
//
//
//        cout << "ghostcounter is " << ghostCounter <<"\n";
        m_num_edges=m_num_edges+ExtraEdges;
    }

    if(Testcase ==45){   //PARTIAL CURVED DAM BREAK
        int ExtraEdges = 40;
        int ghostCounter = 0;
        imatrix EdgeInfo_TMP(7,m_num_edges);

        for(int is=1;is<=m_num_edges;++is){
                for (int i=1;i<=7;i++){
                    EdgeInfo_TMP(i,is) = EdgeInfo(i,is);
                }
        }

        EdgeInfo.resize(7,m_num_edges+ExtraEdges);

        for(int is=1;is<=m_num_edges;++is){
                for (int i=1;i<=7;i++){
                      EdgeInfo(i,is) = EdgeInfo_TMP(i,is);
                }
        }
//
//
        for(int is=1;is<=m_num_edges;++is){
            if( (EdgeInfo(3,is) % 20 == 2) && (EdgeInfo(3,is) % 40 != 2) && (EdgeInfo(4,is) % 20 == 3)&& (EdgeInfo(4,is) % 40 != 3)){

                ghostCounter = ghostCounter+1;
                EdgeInfo(1,m_num_edges+ghostCounter)= EdgeInfo(1,is);
                EdgeInfo(2,m_num_edges+ghostCounter)= EdgeInfo(2,is);
                EdgeInfo(3,m_num_edges+ghostCounter)= EdgeInfo(4,is);
                EdgeInfo(4,m_num_edges+ghostCounter)= -1;
                EdgeInfo(5,m_num_edges+ghostCounter)= EdgeInfo(6,is);
                EdgeInfo(6,m_num_edges+ghostCounter)= -1;
                EdgeInfo(7,m_num_edges+ghostCounter)= EdgeInfo(7,is);

                EdgeInfo(4,is) = -1 ;
                EdgeInfo(6,is)=-1;

            }
        }
//        cout << "ghostcounter is " << ghostCounter <<"\n";
        ghostCounter = 1;
        for(int ie=24;ie<=1584;ie=ie+40){

            ElementToEdge(4,ie)   = m_num_edges + ghostCounter;
            ElemEdgeMasterSlave(4,ie)=+1;
            ElemEdgeOrientation(4,ie)=1;   //order is never reversed for left element! as it is master

            ghostCounter = ghostCounter+1;

        }
//
//
//        cout << "ghostcounter is " << ghostCounter <<"\n";
        m_num_edges=m_num_edges+ExtraEdges;
    }









    if(Testcase ==20){   //PARTIAL DAM BREAK, VARIABLE NelemX,NelemY
        int ExtraEdges= int(0.9*NelemX); // the dam is supposed to have length of 1, so 1/10 of the elements
        int CutEdges = int(0.1*NelemX);
        int ghostCounter = 0;
        imatrix EdgeInfo_TMP(7,m_num_edges);

        for(int is=1;is<=m_num_edges;++is){
                for (int i=1;i<=7;i++){
                    EdgeInfo_TMP(i,is) = EdgeInfo(i,is);
                }
        }

        EdgeInfo.resize(7,m_num_edges+ExtraEdges);

        for(int is=1;is<=m_num_edges;++is){
                for (int i=1;i<=7;i++){
                      EdgeInfo(i,is) = EdgeInfo_TMP(i,is);
                }
        }
//
//
//        for(int is=1;is<=m_num_edges;++is){

         int startIndex = (NelemX+1)*NelemY + NelemY/2 +1;
         int endIndex = m_num_edges - NelemY/2 ;
         int increment =NelemX+1;

         int MidIndex = startIndex + increment*(ExtraEdges/2-1);
         int startIndex2 = MidIndex + (CutEdges+1) * increment;

//
//        cout << " startIndex: " << startIndex <<"\n";
//        cout << " endIndex: " << endIndex<<"\n";
//        cout << " increment: " << increment<<"\n";
//        cout << " MidIndex: " << MidIndex<<"\n";
//        cout << " startIndex2: " << startIndex2<<"\n";
//        for(int is=1661;is<=3260;is+=41){       //40x40 case
        for(int is=startIndex;is<=MidIndex;is+=increment){
                ghostCounter = ghostCounter+1;
                EdgeInfo(1,m_num_edges+ghostCounter)= EdgeInfo(1,is);
                EdgeInfo(2,m_num_edges+ghostCounter)= EdgeInfo(2,is);
                EdgeInfo(3,m_num_edges+ghostCounter)= EdgeInfo(4,is);
                EdgeInfo(4,m_num_edges+ghostCounter)= -1;
                EdgeInfo(5,m_num_edges+ghostCounter)= EdgeInfo(6,is);
                EdgeInfo(6,m_num_edges+ghostCounter)= -1;
                EdgeInfo(7,m_num_edges+ghostCounter)= EdgeInfo(7,is);
                EdgeInfo(4,is) = -1 ;
                EdgeInfo(6,is)=-1;
}

        for(int is=startIndex2;is<=endIndex;is+=increment){
                ghostCounter = ghostCounter+1;
                EdgeInfo(1,m_num_edges+ghostCounter)= EdgeInfo(1,is);
                EdgeInfo(2,m_num_edges+ghostCounter)= EdgeInfo(2,is);
                EdgeInfo(3,m_num_edges+ghostCounter)= EdgeInfo(4,is);
                EdgeInfo(4,m_num_edges+ghostCounter)= -1;
                EdgeInfo(5,m_num_edges+ghostCounter)= EdgeInfo(6,is);
                EdgeInfo(6,m_num_edges+ghostCounter)= -1;
                EdgeInfo(7,m_num_edges+ghostCounter)= EdgeInfo(7,is);
                EdgeInfo(4,is) = -1 ;
                EdgeInfo(6,is)=-1;
}
////        cout << "ghostcounter is " << ghostCounter <<"\n";
        ghostCounter = 1;
//                for(int ie=50;ie<=9950;ie=ie+100){
//        for(int ie=21;ie<=1581;ie=ie+40){

            int eleStartIndex = NelemX/2+1;
            int eleEndIndex = m_num_elements - NelemY/2 +1;
            int eleIncrement = NelemX;
            int eleMidIndex = eleStartIndex + eleIncrement*(ExtraEdges/2-1);
            int eleStartIndex2 = eleMidIndex + (CutEdges+1) * eleIncrement;



//        cout << " eleStartIndex: " << eleStartIndex <<"\n";
//        cout << " eleEndIndex: " << eleEndIndex<<"\n";
//        cout << " eleIncrement: " << eleIncrement<<"\n";
//        cout << " eleMidIndex: " << eleMidIndex<<"\n";
//        cout << " eleStartIndex2: " << eleStartIndex2<<"\n";


          for(int ie=eleStartIndex;ie<=eleMidIndex;ie=ie+eleIncrement){
                ElementToEdge(4,ie)   = m_num_edges + ghostCounter;
                ElemEdgeMasterSlave(4,ie)=+1;
                ElemEdgeOrientation(4,ie)=1;   //order is never reversed for left element! as it is master
                ghostCounter = ghostCounter+1;
        }
          for(int ie=eleStartIndex2;ie<=eleEndIndex;ie=ie+eleIncrement){
                ElementToEdge(4,ie)   = m_num_edges + ghostCounter;
                ElemEdgeMasterSlave(4,ie)=+1;
                ElemEdgeOrientation(4,ie)=1;   //order is never reversed for left element! as it is master
                ghostCounter = ghostCounter+1;
        }



//
//
//        cout << "ghostcounter is " << ghostCounter <<"\n";
        m_num_edges=m_num_edges+ExtraEdges;
//        cout << " WE NOW HAVE " << m_num_edges << " total edges!\n";
    }



  NormalsX.resize(ngl*m_num_edges,1);
  NormalsY.resize(ngl*m_num_edges,1);
  Scal.resize(ngl*m_num_edges,1);

  for (int is = 1; is<=m_num_edges;is++){
    for (int i =1; i<=ngl;i++){
        int id = (is-1)*ngl + i;
        NormalsX(id,1)	=   nx_global(EdgeInfo(5,is)*ngl+i,EdgeInfo(3,is)+1);
        NormalsY(id,1)	=   ny_global(EdgeInfo(5,is)*ngl+i,EdgeInfo(3,is)+1);
        Scal(id,1) 	    =   scal_global(EdgeInfo(5,is)*ngl+i,EdgeInfo(3,is)+1);
    }

  }

}






void Mesh::GenerateMesh(const dfloat xL,const dfloat xR,const dfloat yL,const dfloat yR,const bool PeriodicBD_X,const bool PeriodicBD_Y)
{

      int ngl2=ngl*ngl;

    cout <<"NelemX: "<<NelemX<< "  NelemY: "<<NelemY <<".\n";
    m_num_elements=NelemX*NelemY;



    m_num_edges=(NelemX+1)*NelemY + (NelemY+1)*NelemX;

    if (PeriodicBD_X){
        m_num_edges = m_num_edges -NelemX;
    }
    if (PeriodicBD_Y){
        m_num_edges = m_num_edges -NelemY;

    }
    dfloat deltaX=(xR-xL)/NelemX;
    dfloat deltaY=(yR-yL)/NelemY;


    dfloat nx,ny;

    x_global.resize(m_num_elements*ngl*ngl,1);
    y_global.resize(m_num_elements*ngl*ngl,1);
    xXi_global.resize(m_num_elements*ngl*ngl,1);
    xEta_global.resize(m_num_elements*ngl*ngl,1);
    yXi_global.resize(m_num_elements*ngl*ngl,1);
    yEta_global.resize(m_num_elements*ngl*ngl,1);
    J_global.resize(m_num_elements*ngl*ngl,1);
    nx_global.resize(4*ngl,m_num_elements);
    ny_global.resize(4*ngl,m_num_elements);
    scal_global.resize(4*ngl,m_num_elements);

    EdgeInfo.resize(7,m_num_edges);
    //store global edge number for element local sides 1-4
    ElementToEdge.resize(4,m_num_elements);


    dfloat x_xi=0.5*deltaX;
    dfloat y_eta=0.5*deltaY;
    dfloat Jac=x_xi*y_eta;

    //generate mesh
    for(int ieY=0;ieY<NelemY;++ieY){
        for(int ieX=0;ieX<NelemX;++ieX){
          for(int j=0;j<ngl;++j){
            for(int i=0;i<ngl;++i){
                int id = (ieY*NelemX+ieX)*ngl2   +j*ngl+i+1;
                x_global(id) = xL + (x_GL(i+1)+1.0)/2.0 * deltaX + ieX* deltaX;
                y_global(id) = yL + (x_GL(j+1)+1.0)/2.0 * deltaY + ieY* deltaY;
//                cout <<"ID: " <<id <<" (x,y) = ("<<x_global(id)<<" , "<<y_global(id)<<")\n";
            }
          }
        }
    }



    for (unsigned ie = 0; ie < m_num_elements; ++ie){
        for(int j=1;j<=ngl;++j){
            for(int i=1;i<=ngl;++i){

                int id = ie*ngl2   +(j-1)*ngl+i;

                xXi_global(id) = x_xi;
                xEta_global(id) = 0.0;
                yXi_global(id) = 0.0;
                yEta_global(id) = y_eta;
                J_global(id) = Jac;
            }
          }


    }

    int jLoopBound = NelemY;

    if (PeriodicBD_X){
     jLoopBound = jLoopBound - 1;
    }
        //eta==const edges first
    for (int j=0;j<=jLoopBound;j++){
        for (int i=0;i<NelemX;i++){
            int is = (j*NelemX+i)+1;

//            cout << "EDGE: " << is <<"\n";
            if (j==0){
                // WE SWAP THE USUAL CONVENTION OF NUMBERING HERE, SUCH THAT THE LEFT ELEMENT ON EACH EDGE IS NOT -1


                    EdgeInfo(3,is) = i;
                    EdgeInfo(5,is) =0;
                    EdgeInfo(4,is) =-1;
                    EdgeInfo(6,is) =-1;


                for(int iloc=1;iloc<=ngl;++iloc){
                    int id = EdgeInfo(5,is)*ngl  +iloc;
                    ny_global(id,EdgeInfo(3,is)+1) = -1.0;
                    nx_global(id,EdgeInfo(3,is)+1) = 0.0;
                    scal_global(id,EdgeInfo(3,is)+1) = x_xi;
//                     cout << "nx_global (" << id << " , " <<EdgeInfo(3,is)+1 <<" ) : " <<nx_global(id,EdgeInfo(3,is)+1)  <<"\n";
                }


                if (PeriodicBD_X){
                    EdgeInfo(4,is) = (NelemY-1)*NelemX+i;
                    EdgeInfo(6,is) =2;

                    for(int iloc=1;iloc<=ngl;++iloc){
                        int id = EdgeInfo(6,is)*ngl  +iloc;
                        ny_global(id,EdgeInfo(4,is)+1) = 1.0;
                        nx_global(id,EdgeInfo(4,is)+1) = 0.0;
                        scal_global(id,EdgeInfo(4,is)+1) = x_xi;
                    }
                }




            }else if(j==NelemY){


                if (!PeriodicBD_X){

                    EdgeInfo(3,is) =(NelemY-1)*NelemX+i;
                    EdgeInfo(5,is) =2;
                    EdgeInfo(4,is) =-1;
                    EdgeInfo(6,is) =-1;
                    for(int iloc=1;iloc<=ngl;++iloc){
                        int id = EdgeInfo(5,is)*ngl  +iloc;
                        ny_global(id,EdgeInfo(3,is)+1) = 1.0;
                        nx_global(id,EdgeInfo(3,is)+1) = 0.0;
                        scal_global(id,EdgeInfo(3,is)+1) = x_xi;
//                        cout << "nx_global (" << id << " , " <<EdgeInfo(3,is)+1 <<" ) : " <<nx_global(id,EdgeInfo(3,is)+1)  <<"\n";
                    }

                }


//
//                if (PeriodicBD_X){
//
//                    EdgeInfo(4,is) = (NelemY-1)*NelemX+i;
//                    EdgeInfo(6,is) =2;
//
////                    for(int iloc=1;iloc<=ngl;++iloc){
////                        int id = EdgeInfo(6,is)*ngl  +iloc;
////                        ny_global(id,EdgeInfo(4,is)+1) = 1.0;
////                        nx_global(id,EdgeInfo(4,is)+1) = 0.0;
////                        scal_global(id,EdgeInfo(4,is)+1) = x_xi;
////                    }
//                }




            }
            else{
                EdgeInfo(3,is) =(j-1)*NelemX+i;
                EdgeInfo(4,is)= j   *NelemX+i;
                EdgeInfo(5,is) =2;
                EdgeInfo(6,is) =0;
                for(int iloc=1;iloc<=ngl;++iloc){
                    int id = EdgeInfo(5,is)*ngl  +iloc;
                    ny_global(id,EdgeInfo(3,is)+1) = 1.0;
                    nx_global(id,EdgeInfo(3,is)+1) = 0.0;
                    scal_global(id,EdgeInfo(3,is)+1) = x_xi;
                }
                for(int iloc=1;iloc<=ngl;++iloc){
                    int id = EdgeInfo(6,is)*ngl  +iloc;
                    ny_global(id,EdgeInfo(4,is)+1) = -1.0;
                    nx_global(id,EdgeInfo(4,is)+1) = 0.0;
                    scal_global(id,EdgeInfo(4,is)+1) = x_xi;
                }

            };
            EdgeInfo(1,is)=0;
            EdgeInfo(2,is)=0;

            EdgeInfo(7,is) =1;


            ElementToEdge(EdgeInfo(5,is)+1,EdgeInfo(3,is)+1)=is;


            if ( EdgeInfo(4,is) != -1){
//                cout <<"Element: " <<EdgeInfo(4,is)+1 <<"   side 1 is edge  "<<is <<"\n" ;
                ElementToEdge(EdgeInfo(6,is)+1,EdgeInfo(4,is)+1)=is;

            }

        }
    }

    int iLoopBound = NelemX;

    if (PeriodicBD_Y){
     iLoopBound = iLoopBound - 1;
    }
    //xi==const edges now
    for (int j=0;j<NelemY;j++){
        for (int i=0;i<=iLoopBound;i++){
            int is = ((jLoopBound+1)*NelemX) + (j*(iLoopBound+1)+i)+1;        //12 + ....       2 +2
//            cout << "EDGE xi: " << is <<"\n";

            if (i==0){
                // WE SWAP THE USUAL CONVENTION OF NUMBERING HERE, SUCH THAT THE LEFT ELEMENT ON EACH EDGE IS NOT -1



                    EdgeInfo(3,is)  = j*NelemX;
                    EdgeInfo(5,is) =3;
                    EdgeInfo(4,is)  = -1;
                    EdgeInfo(6,is) =-1;

                for(int iloc=1;iloc<=ngl;++iloc){
                    int id = EdgeInfo(5,is)*ngl  +iloc;
                    ny_global(id,EdgeInfo(3,is)+1) = 0.0;
                    nx_global(id,EdgeInfo(3,is)+1) = -1.0;
                    scal_global(id,EdgeInfo(3,is)+1) = y_eta;
                }

                if (PeriodicBD_Y){

                    EdgeInfo(4,is)  =  (j+1)*NelemX-1;
                    EdgeInfo(6,is) =1;
                    for(int iloc=1;iloc<=ngl;++iloc){
                        int id = EdgeInfo(6,is)*ngl  +iloc;
                        ny_global(id,EdgeInfo(4,is)+1) = 0.0;
                        nx_global(id,EdgeInfo(4,is)+1) = 1.0;
                        scal_global(id,EdgeInfo(4,is)+1) = y_eta;
                    }
                }




            }else if(i==NelemX){


            // these edges dont exist if Periodic Y boundaries are chosen
                 if (!PeriodicBD_Y){
                    EdgeInfo(3,is)  =  (j+1)*NelemX-1;
                    EdgeInfo(5,is) =1;
                    EdgeInfo(4,is)  = -1;
                    EdgeInfo(6,is) =-1;
                    for(int iloc=1;iloc<=ngl;++iloc){
                        int id = EdgeInfo(5,is)*ngl  +iloc;
                        ny_global(id,EdgeInfo(3,is)+1) = 0.0;
                        nx_global(id,EdgeInfo(3,is)+1) = 1.0;
                        scal_global(id,EdgeInfo(3,is)+1) = y_eta;
                    }
                 }







            }
            else{
                EdgeInfo(3,is)  =  j*NelemX+i-1;
                EdgeInfo(4,is)  =  j*NelemX+i;
                EdgeInfo(5,is) =1;
                EdgeInfo(6,is) =3;
                for(int iloc=1;iloc<=ngl;++iloc){
                     int id = EdgeInfo(5,is)*ngl  +iloc;
                    ny_global(id,EdgeInfo(3,is)+1) = 0.0;
                    nx_global(id,EdgeInfo(3,is)+1) = 1.0;
                    scal_global(id,EdgeInfo(3,is)+1) = y_eta;
                }
                for(int iloc=1;iloc<=ngl;++iloc){
                    int id = EdgeInfo(6,is)*ngl  +iloc;
                    ny_global(id,EdgeInfo(4,is)+1) = 0.0;
                    nx_global(id,EdgeInfo(4,is)+1) = -1.0;
                    scal_global(id,EdgeInfo(4,is)+1) = y_eta;
                }


            };
            EdgeInfo(1,is)=0;
            EdgeInfo(2,is)=0;

            EdgeInfo(7,is) =1;

//            cout <<"Element: " <<EdgeInfo(3,is)+1 <<"   side 2 is edge  "<<is <<"\n" ;
                ElementToEdge(EdgeInfo(5,is)+1,EdgeInfo(3,is)+1)=is;


            if ( EdgeInfo(4,is) != -1){
//            cout <<"Element: " <<EdgeInfo(4,is)+1 <<"   side 4 is edge  "<<is <<"\n" ;
                ElementToEdge(EdgeInfo(6,is)+1,EdgeInfo(4,is)+1)=is;

            }
        }
    }

//    for (unsigned ie = 1; ie <= m_num_elements; ++ie){
//        cout <<"Element: "<<ie<<"\n";
//        for (unsigned is=1;is<=4;is++){
//            cout <<"side "<<is<<" is global edge: "<<ElementToEdge(is,ie)<<"\n";
//            cout <<"global edge "<< ElementToEdge(is,ie) <<" has left element "<< EdgeInfo(3,ElementToEdge(is,ie)) <<" and right element "<< EdgeInfo(4,ElementToEdge(is,ie)) <<"\n";
//        }
//
//    }

//for (int is=1;is<=m_num_edges;is++){
//
//    cout << "Edge Number: "<<is<<"   ";
//    cout <<"Left Element: " <<EdgeInfo(3,is)+1 <<"  Right Element:  "<<EdgeInfo(4,is)+1 <<"    " ;
//    cout <<"Edge Left: " <<EdgeInfo(5,is) <<"  Edge Right:  "<<EdgeInfo(6,is)<<"\n" ;
//    cout <<"nx: " ;
//    for(int iloc=1;iloc<=ngl;++iloc){
//        int id = (EdgeInfo(5,is))*ngl  +iloc;
//        cout << " "<<nx_global(id,EdgeInfo(3,is)+1);
////        cout << "nx_global (" << id << " , " <<EdgeInfo(3,is)+1 <<" ) : " <<nx_global(id,EdgeInfo(3,is)+1)  <<"\n";
//    }
//    cout <<"      ny: " ;
//    for(int iloc=1;iloc<=ngl;++iloc){
//        int id = (EdgeInfo(5,is))*ngl  +iloc;
//        cout << " "<<ny_global(id,EdgeInfo(3,is)+1);
//    }
//    cout <<"\n";
//}



        // STORING EDGE INFORMATION!!
    //1 = start node id of edge
    //2 = end node id of edge
    //3 = id of left element
    //4 = id of right element
    //5 = local edge left element
    //6 = local edge right element
    //7 = is orientation swapped?

}

void Mesh::ReadMesh(const string meshFile)
{

//,dfloat T, dfloat g_const

    corners.resize(4,2);
    x_phy.resize(ngl,ngl);
    y_phy.resize(ngl,ngl);
    x_xi.resize(ngl,ngl);
    x_eta.resize(ngl,ngl);
    y_xi.resize(ngl,ngl);
    y_eta.resize(ngl,ngl);
    J.resize(ngl,ngl);
    x_bndy.resize(ngl,4);
    y_bndy.resize(ngl,4);
    scal.resize(ngl,4);
    nx.resize(ngl,4);
    ny.resize(ngl,4);




    std::ifstream InputStream;
    string filename=meshFile;
    InputStream.open(filename.c_str());

    if (!InputStream)
    {
      std::string error_message("ERROR: Mesh file not found: ");
      error_message += filename;
      throw std::invalid_argument(error_message);
    }

    std::string FileFormat;
    std::string current_string;
    std::getline(InputStream, current_string);
    std::stringstream current_line(current_string);


    current_line>>FileFormat;
    cout << "Fileformat: "<<FileFormat<<"\n";

    std::getline(InputStream, current_string);
    current_line.clear();
    current_line.str(current_string);


    if (!(current_line >> m_num_nodes >> m_num_edges >> m_num_elements >> m_order_of_boundary_edges))
    {
      std::string error_message("ERROR: Header incorrect! ");
      error_message += filename;
      throw std::invalid_argument(error_message);
    }

    dfloat x_nodes[m_num_nodes];
    dfloat y_nodes[m_num_nodes];
    //initialise chebyshev nodes immediately
    x_cheby.resize(m_order_of_boundary_edges+1);
    w_bary.resize(m_order_of_boundary_edges+1);
    for (unsigned k = 0; k <= m_order_of_boundary_edges; ++k)
    {x_cheby(k+1) = -cos(k * PI / m_order_of_boundary_edges);}  // k = m_order_of_boundary_edges/2 => x_cheby(k+1) = -cos(PI/2) = 0
    //calculate barycentric weights for chebyshev nodes
    BarycentricWeights();

//    for (unsigned k = 1; k <= m_order_of_boundary_edges+1; ++k){
//        cout << "Bary Weights ("<<k<<") : "<<w_bary(k) <<"\n";
//
//   }
//    for (unsigned k = 1; k <= m_order_of_boundary_edges+1; ++k){
//        cout << "x_cheby ("<<k<<") : "<<x_cheby(k) <<"\n";
//
//   }



    for (unsigned i = 0; i < m_num_nodes; ++i)
    {
      std::getline(InputStream, current_string);
      current_line.clear();
      current_line.str(current_string);
      if (!(current_line >> x_nodes[i] >> y_nodes[i]))
           {
              std::string error_message("ERROR: Cant read in Nodes! ");
              error_message += filename;
              throw std::invalid_argument(error_message);
        }
    }


    // STORING EDGE INFORMATION!!
    //1 = start node id of edge
    //2 = end node id of edge
    //3 = id of left element
    //4 = id of right element
    //5 = local edge left element
    //6 = local edge right element
    //7 = is orientation swapped?
//    int EdgeData[7*m_num_edges];

    EdgeInfo.resize(7,m_num_edges);
    //store global edge number for element local sides 1-4
    ElementToEdge.resize(4,m_num_elements);

    for (int is = 1; is <= m_num_edges; ++is)
    {
      std::getline(InputStream, current_string);
      current_line.clear();
      current_line.str(current_string);
      int start_node_id, end_node_id, element_id_on_left, element_id_on_right, side_of_left_element;
      int side_of_right_element;
      if (!(current_line >> start_node_id >> end_node_id
                         >> element_id_on_left >> element_id_on_right
                         >> side_of_left_element >> side_of_right_element))
      {
              std::string error_message("ERROR: Cant read in Edge Stuff! ");
              error_message += filename;
              throw std::invalid_argument(error_message);
      }
      EdgeInfo(1,is)=start_node_id-1;
      EdgeInfo(2,is)= end_node_id-1;
      EdgeInfo(3,is)= element_id_on_left-1;
      EdgeInfo(4,is)= element_id_on_right-1;
      EdgeInfo(5,is)= side_of_left_element-1;

      ElementToEdge(side_of_left_element,element_id_on_left)=is;
//      ElementToEdge(EdgeInfo(5,is)+1,EdgeInfo(3,is)+1)=is;




//      EdgeData[7*is+0]= start_node_id-1;
//      EdgeData[7*is+1]= end_node_id-1;
//      EdgeData[7*is+2]= element_id_on_left-1;
//      EdgeData[7*is+3]= element_id_on_right-1;
//      EdgeData[7*is+4]= side_of_left_element-1;


      // sgn(side_of_right_element)

        int index_direction;
      if (side_of_right_element<0){
            //SIDE IS REVERSED
            index_direction=0;
            EdgeInfo(6,is)= -side_of_right_element-1;
        }else{
            index_direction=1;
            EdgeInfo(6,is)= side_of_right_element-1;
        }
        EdgeInfo(7,is)= index_direction;

      if (EdgeInfo(4,is)!=-1){
          if (side_of_right_element<0){
                //SIDE IS REVERSED
                ElementToEdge(-side_of_right_element,element_id_on_right)=is;
            }else{

                ElementToEdge(side_of_right_element,element_id_on_right)=is;
            }
        }


//          ElementToEdge(EdgeInfo(6,is)+1,EdgeInfo(4,is)+1)=is;


    }



//    for (unsigned ie = 1; ie <= m_num_elements; ++ie){
//        cout <<"Element: "<<ie<<"\n";
//        for (unsigned is=1;is<=4;is++){
//            cout <<"side "<<is<<" is global edge: "<<ElementToEdge(is,ie)<<"\n";
//            cout <<"global edge "<< ElementToEdge(is,ie) <<" has left element "<< EdgeInfo(3,ElementToEdge(is,ie)) <<" and right element "<< EdgeInfo(4,ElementToEdge(is,ie)) <<"\n";
//        }
//
//    }


    // STORING ELEMENT DATA INFORMATION
    //1 = Corner node 1     ID!!
    //2 = Corner node 2     ID!!
    //3 = Corner node 3     ID!!
    //4 = Corner node 4     ID!!

    // STORING ELEMENT Boundary Name INFORMATION
    //1 = Boundary Name side 1
    //2 = Boundary Name side 2
    //3 = Boundary Name side 3
    //4 = Boundary Name side 4

//    unsigned eleInfoNo=4;
//    int ElementData[eleInfoNo*m_num_elements];

    string ElementBoundaryNames[4*m_num_elements];
    Gamma1.resize(m_order_of_boundary_edges+1,2);
    Gamma2.resize(m_order_of_boundary_edges+1,2);
    Gamma3.resize(m_order_of_boundary_edges+1,2);
    Gamma4.resize(m_order_of_boundary_edges+1,2);

    x_global.resize(m_num_elements*ngl*ngl,1);
    y_global.resize(m_num_elements*ngl*ngl,1);
    xXi_global.resize(m_num_elements*ngl*ngl,1);
    xEta_global.resize(m_num_elements*ngl*ngl,1);
    yXi_global.resize(m_num_elements*ngl*ngl,1);
    yEta_global.resize(m_num_elements*ngl*ngl,1);
    J_global.resize(m_num_elements*ngl*ngl,1);
    nx_global.resize(4*ngl,m_num_elements);
    ny_global.resize(4*ngl,m_num_elements);
    scal_global.resize(4*ngl,m_num_elements);


//    double GammaCurvesX[m_num_elements*4*(m_order_of_boundary_edges+1)];
//    double GammaCurvesY[m_num_elements*4*(m_order_of_boundary_edges+1)];
//    for (unsigned i=0;i<m_num_elements*4*(m_order_of_boundary_edges+1);i++){
//        GammaCurvesX[i]=0.0;
//        GammaCurvesY[i]=0.0;





    for (unsigned ie = 0; ie < m_num_elements; ++ie)
    {
       bool curved=false;
       for (unsigned i=1;i<(m_order_of_boundary_edges+1);i++){
            Gamma1(i,1)=0.0;
            Gamma1(i,2)=0.0;
            Gamma2(i,1)=0.0;
            Gamma2(i,2)=0.0;
            Gamma3(i,1)=0.0;
            Gamma3(i,2)=0.0;
            Gamma4(i,1)=0.0;
            Gamma4(i,2)=0.0;
        }
        int CornerIDs[4];
      // Get the corner nodes.
      std::getline(InputStream, current_string);
      current_line.clear();
      current_line.str(current_string);
      if (!(current_line >> CornerIDs[0 ]
                         >> CornerIDs[1 ]
                         >> CornerIDs[2 ]
                         >> CornerIDs[3 ]))

      {
        std::string error_message("ERROR: Cant read in Element Corners! ");
              error_message += filename;
              throw std::invalid_argument(error_message);
      }


      // ISM2 format starts counting from one, so we substract one since we count from zero.
      for (unsigned j = 0; j < 4; ++j){
        CornerIDs[j ]-=1;
//        ElementData[eleInfoNo*ie + j ] -= 1;
      }

      // Get the information about whether the boundary curves are given.
      std::getline(InputStream, current_string);
      current_line.clear();
      current_line.str(current_string);
      int is_gamma_given[4];
      if (!(current_line >> is_gamma_given[0] >> is_gamma_given[1] >> is_gamma_given[2] >> is_gamma_given[3]))
          {
            std::string error_message("ERROR: Cant read in Gamma information! ");
                  error_message += filename;
                  throw std::invalid_argument(error_message);
          }

      // If none of the gammas are given, then we have a fully straight sided quadrilateral.
      // We do not need to set the gammas in this case.
      if (!(is_gamma_given[0] || is_gamma_given[1] || is_gamma_given[2] || is_gamma_given[3])){
        curved=false;}//ElementData[eleInfoNo*i + 4 ]=1;    // side is straight sided
      else
      {
        curved=true;}//ElementData[eleInfoNo*i + 4 ]=0;    // side is curved

        // Determine if the boundary curve is given.
        // If it is not, then we need to build the boundary curve from the two nodes.
        // If it is, simply copy over the boundary curve information to the intermediate storage class.
        for (unsigned j = 0; j < 4; ++j)
        {
            //initialise the corners fmatrix or set it to the current corner nodes
            corners(j+1,1)=    x_nodes[CornerIDs[j]];
            corners(j+1,2)=    y_nodes[CornerIDs[j]];


          //if we dont have a given Gamma curve, we get one by interpolation (chebychev)
          if (is_gamma_given[j] == false)
          {

              unsigned corner_initial;
              unsigned corner_final;

              switch(j){
              case 0:
                corner_initial= CornerIDs[ 0 ];
                corner_final= CornerIDs[ 1 ];
               break;
                case 1:
                corner_initial= CornerIDs[ 1 ];
                corner_final= CornerIDs[ 2 ];
                 break;
              case 2:
                corner_initial= CornerIDs[ 3 ];
                corner_final= CornerIDs[ 2 ];
                 break;
              case 3:
                corner_initial= CornerIDs[ 0 ];
                corner_final= CornerIDs[ 3 ];
                 break;

              }


              double x_initial = x_nodes[corner_initial];
              double x_final = x_nodes[corner_final];

              double y_initial = y_nodes[corner_initial];
              double y_final = y_nodes[corner_final];



            for (unsigned k = 0; k <= m_order_of_boundary_edges; ++k)
            {


              double interpolated_x = ((x_final - x_initial) * x_cheby(k+1) + x_initial + x_final)/2.0;
              double interpolated_y = ((y_final - y_initial) * x_cheby(k+1) + y_initial + y_final)/2.0;


              switch(j){
              case 0:
                Gamma1(k+1,1)= interpolated_x;
                Gamma1(k+1,2)= interpolated_y;
               break;
                case 1:
                Gamma2(k+1,1)= interpolated_x;
                Gamma2(k+1,2)= interpolated_y;
                 break;
              case 2:
                Gamma3(k+1,1)= interpolated_x;
                Gamma3(k+1,2)= interpolated_y;
                 break;
              case 3:
                Gamma4(k+1,1)= interpolated_x;
                Gamma4(k+1,2)= interpolated_y;
                 break;

              }
            }


          }
          else // is_gamma_given is true
          {
            for (unsigned k = 0; k <= m_order_of_boundary_edges; ++k)
            {
              std::getline(InputStream, current_string);
              current_line.clear();
              current_line.str(current_string);
//              if (!(current_line >> GammaCurvesX[i*4*(m_order_of_boundary_edges+1)+j*(m_order_of_boundary_edges+1)+k]
//                                 >> GammaCurvesY[i*4*(m_order_of_boundary_edges+1)+j*(m_order_of_boundary_edges+1)+k]))
              switch(j){
              case 0:
               if (!(current_line >> Gamma1(k+1,1)>> Gamma1(k+1,2)))          {
                std::string error_message("ERROR: Cant read in given Gamma1 Curve! ");
                  error_message += filename;
                  cout << "Error in Line : "<<k <<"\n";
                  throw std::invalid_argument(error_message);
                }
               break;
                case 1:
                if (!(current_line >> Gamma2(k+1,1)>> Gamma2(k+1,2)))          {
                std::string error_message("ERROR: Cant read in given Gamma2 Curve! ");
                  error_message += filename;
                  cout << "Error in Line : "<<k <<"\n";
                  throw std::invalid_argument(error_message);
                }
                 break;
              case 2:
                if (!(current_line >> Gamma3(k+1,1)>> Gamma3(k+1,2)))          {
                std::string error_message("ERROR: Cant read in given Gamma3 Curve! ");
                  error_message += filename;
                  cout << "Error in Line : "<<k <<"\n";
                  throw std::invalid_argument(error_message);
                }
                 break;
              case 3:
                if (!(current_line >> Gamma4(k+1,1)>> Gamma4(k+1,2)))          {
                std::string error_message("ERROR: Cant read in given Gamma4 Curve! ");
                  error_message += filename;
                  cout << "Error in Line : "<<k <<"\n";
                  throw std::invalid_argument(error_message);
                }
                 break;

              }


            }
          }
        }


      // Extract the boundary names
      std::getline(InputStream, current_string);
      current_line.clear();
      current_line.str(current_string);
      for (unsigned j = 0; j < 4; ++j){
        if (!(current_line >> ElementBoundaryNames[4*ie +j]))
          {
                std::string error_message("ERROR: Cant read in given Boundary Names! ");
                  error_message += filename;
                  throw std::invalid_argument(error_message);
          }}



      ConstructMappedGeometry(curved);


      for(int j=1;j<=ngl;++j){
        for(int i=1;i<=ngl;++i){

            int id = ie*ngl*ngl   +(j-1)*ngl+i;
            x_global(id) = x_phy(i,j);
            y_global(id) = y_phy(i,j);
            xXi_global(id) = x_xi(i,j);
            xEta_global(id) = x_eta(i,j);
            yXi_global(id) = y_xi(i,j);
            yEta_global(id) = y_eta(i,j);
            J_global(id) = J(i,j);
        }
      }

      for(int is=1;is<=4;++is){
        for(int i=1;i<=ngl;++i){

            int id = (is-1)*ngl  +i;

            nx_global(id,ie+1) = nx(i,is);
            ny_global(id,ie+1) = ny(i,is);
            scal_global(id,ie+1) = scal(i,is);
        }
      }



    }//ele loop




}






void Mesh :: TransfiniteQuadMetrics(const dfloat psi,const dfloat eta,dfloat *X_psi,dfloat *X_eta,dfloat *Y_psi,dfloat *Y_eta){
// Computation of the metric terms on a curve-bounded quadrilateral


      fmatrix Xref,Xcomp,Xpcomp;
      Xref.resize(4,2);
      Xcomp.resize(4,2);
      Xpcomp.resize(4,2);

      EvaluateAt(Gamma1,-1.0,&Xref(1,1),&Xref(1,2));
      EvaluateAt(Gamma1, 1.0,&Xref(2,1),&Xref(2,2));
      EvaluateAt(Gamma3, 1.0,&Xref(3,1),&Xref(3,2));
      EvaluateAt(Gamma3,-1.0,&Xref(4,1),&Xref(4,2));

      EvaluateAt(Gamma1,psi,&Xcomp(1,1),&Xcomp(1,2));
      EvaluateAt(Gamma2,eta,&Xcomp(2,1),&Xcomp(2,2));
      EvaluateAt(Gamma3,psi,&Xcomp(3,1),&Xcomp(3,2));
      EvaluateAt(Gamma4,eta,&Xcomp(4,1),&Xcomp(4,2));

      DerivativeAt(Gamma1,psi,&Xpcomp(1,1),&Xpcomp(1,2));
      DerivativeAt(Gamma2,eta,&Xpcomp(2,1),&Xpcomp(2,2));
      DerivativeAt(Gamma3,psi,&Xpcomp(3,1),&Xpcomp(3,2));
      DerivativeAt(Gamma4,eta,&Xpcomp(4,1),&Xpcomp(4,2));

      *X_psi = 0.5*(Xcomp(2,1)-Xcomp(4,1)+(1.0-eta)*Xpcomp(1,1)+(1.0+eta)*Xpcomp(3,1)) -  0.25*((1.0-eta)*(Xref(2,1)-Xref(1,1))+(1.0+eta)*(Xref(3,1)-Xref(4,1)));

      *Y_psi = 0.5*(Xcomp(2,2)-Xcomp(4,2)+(1.0-eta)*Xpcomp(1,2)+(1.0+eta)*Xpcomp(3,2)) -0.25*((1.0-eta)*(Xref(2,2)-Xref(1,2))+(1.0+eta)*(Xref(3,2)-Xref(4,2)));

      *X_eta = 0.5*((1.0-psi)*Xpcomp(4,1)+(1.0+psi)*Xpcomp(2,1)+Xcomp(3,1)-Xcomp(1,1)) -  0.25*((1.0-psi)*(Xref(4,1)-Xref(1,1))+(1.0+psi)*(Xref(3,1)-Xref(2,1)));

      *Y_eta = 0.5*((1.0-psi)*Xpcomp(4,2)+(1.0+psi)*Xpcomp(2,2)+Xcomp(3,2)-Xcomp(1,2)) - 0.25*((1.0-psi)*(Xref(4,2)-Xref(1,2))+(1.0+psi)*(Xref(3,2)-Xref(2,2)));

}



void Mesh :: TransfiniteQuadMap(const dfloat psi,const dfloat eta,dfloat *x_out,dfloat *y_out){
// Mapping of the reference square to a curve-bounded quadrilateral
//      IMPLICIT NONE
//      TYPE(CurveInterpolant),DIMENSION(4),INTENT(IN)  :: GammaCurves
//      REAL(KIND=RP)                      ,INTENT(IN)  :: psi,eta
//      REAL(KIND=RP)                      ,INTENT(OUT) :: x_out,y_out

//      REAL(KIND=RP)         ,DIMENSION(4,2)           :: Xref,Xcomp
      dfloat x1_ref, y1_ref,x1_cp, y1_cp;
      dfloat x2_ref, y2_ref,x2_cp, y2_cp;
      dfloat x3_ref, y3_ref,x3_cp, y3_cp;
      dfloat x4_ref, y4_ref,x4_cp, y4_cp;

      EvaluateAt(Gamma1,-1.0,&x1_ref, &y1_ref);
      EvaluateAt(Gamma1, 1.0,&x2_ref, &y2_ref);
      EvaluateAt(Gamma3, 1.0,&x3_ref, &y3_ref);
      EvaluateAt(Gamma3,-1.0,&x4_ref, &y4_ref);

      EvaluateAt(Gamma1,psi,&x1_cp, &y1_cp);
      EvaluateAt(Gamma2,eta,&x2_cp, &y2_cp);
      EvaluateAt(Gamma3,psi,&x3_cp, &y3_cp);
      EvaluateAt(Gamma4,eta,&x4_cp, &y4_cp);

      *x_out = 1.0/2.0*((1.0-psi)*x4_cp+(1.0+psi)*x2_cp+(1.0-eta)*x1_cp+(1.0+eta)*x3_cp) - 1.0/4.0*((1.0-psi)*((1.0-eta)*x1_ref+(1.0+eta)*x4_ref)+(1.0+psi)*((1.0-eta)*x2_ref+(1.0+eta)*x3_ref));

      *y_out = 1.0/2.0*((1.0-psi)*y4_cp+(1.0+psi)*y2_cp+(1.0-eta)*y1_cp+(1.0+eta)*y3_cp) - 1.0/4.0*((1.0-psi)*((1.0-eta)*y1_ref+(1.0+eta)*y4_ref)  +(1.0+psi)*((1.0-eta)*y2_ref+(1.0+eta)*y3_ref));

}

void Mesh :: QuadMap(const dfloat psi,const dfloat eta,dfloat *x_out,dfloat *y_out){
// Mapping of the reference square to a straight sided quadrilateral

      *x_out = 0.25*(corners(1,1)*(1.0-psi)*(1.0-eta)+corners(2,1)*(1.0+psi)*(1.0-eta)+ corners(3,1)*(1.0+psi)*(1.0+eta)+corners(4,1)*(1.0-psi)*(1.0+eta));
      *y_out = 0.25*(corners(1,2)*(1.0-psi)*(1.0-eta)+corners(2,2)*(1.0+psi)*(1.0-eta)+ corners(3,2)*(1.0+psi)*(1.0+eta)+corners(4,2)*(1.0-psi)*(1.0+eta));

}







void Mesh :: QuadMapMetrics(const dfloat xi,const dfloat eta,dfloat *X_psi,dfloat *X_eta,dfloat *Y_psi,dfloat *Y_eta){
// Metric terms on a straight sided quadrilateral

//we assume we have the correct corners stored in the corners fmatrix
      *X_psi = 0.25*((1.0-eta)*(corners(2,1)-corners(1,1))+(1.0+eta)*(corners(3,1)-corners(4,1)));
      *Y_psi = 0.25*((1.0-eta)*(corners(2,2)-corners(1,2))+(1.0+eta)*(corners(3,2)-corners(4,2)));
      *X_eta = 0.25*((1.0-xi)*(corners(4,1)-corners(1,1))+(1.0+xi)*(corners(3,1)-corners(2,1)));
      *Y_eta = 0.25*((1.0-xi)*(corners(4,2)-corners(1,2))+(1.0+xi)*(corners(3,2)-corners(2,2)));

}



//const double GammaX1[],const double GammaX2[],const double GammaX3[],const double GammaX4[],
//                                     const double GammaY1[],const double GammaY2[],const double GammaY3[],const double GammaY4[],
void Mesh :: ConstructMappedGeometry(const bool Curved){
// Constructor of geometry and metric terms for quadrilateral domains


dfloat xXi,xEta,yXi,yEta,Jtemp;
dfloat x_tmp, y_tmp;

     for (int j = 1;j<=ngl;j++){
         for (int i = 1;i<=ngl;i++){
            if (Curved){
               TransfiniteQuadMap(x_GL(i),x_GL(j),&x_phy(i,j),&y_phy(i,j));
               TransfiniteQuadMetrics(x_GL(i),x_GL(j),&x_xi(i,j),&x_eta(i,j),&y_xi(i,j),&y_eta(i,j));
            }else{
               QuadMap(x_GL(i),x_GL(j),&x_phy(i,j),&y_phy(i,j));
               QuadMapMetrics(x_GL(i),x_GL(j),&x_xi(i,j),&x_eta(i,j),&y_xi(i,j),&y_eta(i,j));
                }
            J(i,j)= x_xi(i,j)*y_eta(i,j)-x_eta(i,j)*y_xi(i,j);
        }
     }


//     for (int j = 1;j<=ngl;j++){
//         for (int i = 1;i<=ngl;i++){
//            cout <<"(x,y) = (" << x_phy(i,j)<< "," << y_phy(i,j) <<")\n";
//        }
//        cout <<"\n";
//     }


     for (int j = 1;j<=ngl;j++){
            if (Curved){
                TransfiniteQuadMap(1.0,x_GL(j),&x_bndy(j,2),&y_bndy(j,2));
                TransfiniteQuadMetrics(1.0,x_GL(j),&xXi,&xEta,&yXi,&yEta);
            }else{
                QuadMap(1.0,x_GL(j),&x_bndy(j,2),&y_bndy(j,2));
                QuadMapMetrics(1.0,x_GL(j),&xXi,&xEta,&yXi,&yEta);
                }
         Jtemp = xXi*yEta-xEta*yXi;
         scal(j,2)      = sqrt(yEta*yEta + xEta*xEta);
         nx(j,2) = copysign(1.0,Jtemp)*(yEta/scal(j,2));
         ny(j,2) = copysign(1.0,Jtemp)*(-xEta/scal(j,2));
             if (Curved){
                TransfiniteQuadMap(-1.0,x_GL(j),&x_bndy(j,4),&y_bndy(j,4));
                TransfiniteQuadMetrics(-1.0,x_GL(j),&xXi,&xEta,&yXi,&yEta);
            }else{
                QuadMap(-1.0,x_GL(j),&x_bndy(j,4),&y_bndy(j,4));
                QuadMapMetrics(-1.0,x_GL(j),&xXi,&xEta,&yXi,&yEta);
            }
         Jtemp = xXi*yEta-xEta*yXi;
         scal(j,4)      =  sqrt(yEta*yEta + xEta*xEta);
         nx(j,4) = -copysign(1.0,Jtemp)*(yEta/scal(j,4));
         ny(j,4) = -copysign(1.0,Jtemp)*(-xEta/scal(j,4));
     }
     for (int i = 1;i<=ngl;i++){
            if (Curved){
                TransfiniteQuadMap(x_GL(i),-1.0,&x_bndy(i,1),&y_bndy(i,1));
                TransfiniteQuadMetrics(x_GL(i),-1.0,&xXi,&xEta,&yXi,&yEta);
            }else{
                QuadMap(x_GL(i),-1.0,&x_bndy(i,1),&y_bndy(i,1));
                QuadMapMetrics(x_GL(i),-1.0,&xXi,&xEta,&yXi,&yEta);
            }
         Jtemp = xXi*yEta-xEta*yXi;
         scal(i,1)      =  sqrt(yXi*yXi + xXi*xXi);
         nx(i,1) = -copysign(1.0,Jtemp)*(-yXi/scal(i,1));
         ny(i,1) = -copysign(1.0,Jtemp)*(xXi/scal(i,1));
            if (Curved){
                TransfiniteQuadMap(x_GL(i),1.0,&x_bndy(i,3),&y_bndy(i,3));
                TransfiniteQuadMetrics(x_GL(i),1.0,&xXi,&xEta,&yXi,&yEta);
            }else{
                QuadMap(x_GL(i),1.0,&x_bndy(i,3),&y_bndy(i,3));
                QuadMapMetrics(x_GL(i),1.0,&xXi,&xEta,&yXi,&yEta);
            }
         Jtemp = xXi*yEta-xEta*yXi;
         scal(i,3)      = sqrt(yXi*yXi + xXi*xXi);
         nx(i,3) = copysign(1.0,Jtemp)*(-yXi/scal(i,3));
         ny(i,3) = copysign(1.0,Jtemp)*(xXi/scal(i,3));
     }

}























void Mesh :: BarycentricWeights(){


        for(int i=1;i<=m_order_of_boundary_edges+1;i++){
        w_bary(i)=1.0;
        };


        for (int j=2;j<=m_order_of_boundary_edges+1;j++){
             for (int k=1;k<j;k++){
                w_bary(k)=w_bary(k)*(x_cheby(k)-x_cheby(j));
                w_bary(j)=w_bary(j)*(x_cheby(j)-x_cheby(k));
            };
        };

        for (int j=1;j<=m_order_of_boundary_edges+1;j++){
         w_bary(j)=1.0/w_bary(j);
        };

};




void Mesh :: EvaluateAt(const fmatrix Gamma,const dfloat s,dfloat *x_point,dfloat *y_point){
// Evaluates a member of the CurveInterpolant type at a point s
dfloat GammaX[ngl],GammaY[ngl];

for (int i = 0; i<ngl;i++){
    GammaX[i] = Gamma(i+1,1);
    GammaY[i] = Gamma(i+1,2);

}
      LagrangeInterpolation(s,GammaX,x_point);
      LagrangeInterpolation(s,GammaY,y_point);
}

void Mesh :: DerivativeAt(const fmatrix Gamma,const dfloat s,dfloat *x_point_prime,dfloat *y_point_prime){
// Evaluates the derivative of a member of the CurveInterpolant type at a point s
dfloat GammaX[ngl],GammaY[ngl];

for (int i = 0; i<ngl;i++){
    GammaX[i] = Gamma(i+1,1);
    GammaY[i] = Gamma(i+1,2);

}

LagrangeInterpolantDerivative(s,GammaX,x_point_prime);
LagrangeInterpolantDerivative(s,GammaY,y_point_prime);

}



void Mesh :: LagrangeInterpolantDerivative(const dfloat xpt,const dfloat functionvals[],dfloat *p_prime){

      bool atNode      = false;
      dfloat numerator   = 0.0;
      dfloat denominator=0.0;
      dfloat p,t;
      int k;

     for (int j = 0; j<=m_order_of_boundary_edges;j++){
        if(fabs(xpt-x_cheby(j+1))<pow(10.0,-14)){
            atNode = true;
            p = functionvals[j];
            denominator = -1.0*w_bary(j+1);
            k = j;

         }
     }
      if (atNode){
        for (int j = 0; j<=m_order_of_boundary_edges;j++){
            if(j!=k){
               numerator = numerator + w_bary(j+1)*(p-functionvals[j])/(xpt-x_cheby(j+1));
            }
        }
        }else{
             denominator = 0.0;
             LagrangeInterpolation(xpt,functionvals,&p);
             for (int j = 0; j<=m_order_of_boundary_edges;j++){
                t = w_bary(j+1)/(xpt-x_cheby(j+1));
                numerator +=  t*(p-functionvals[j])/(xpt-x_cheby(j+1));
                denominator += t;
             }
         }
      *p_prime = numerator/denominator;

}


void Mesh :: LagrangeInterpolation(const dfloat xpt,const dfloat functionvals[],dfloat *output){
// Barycentric two formulation of Lagrange interpolant



      dfloat numerator,denominator,t;
        bool var1=false;

      numerator   = 0.0;
      denominator = 0.0;

     for (int j = 0; j<=m_order_of_boundary_edges;j++){
         if(fabs(xpt-x_cheby(j+1))<pow(10.0,-14)){
            var1=true;
            *output = functionvals[j];
            break;
         }else{
            t = w_bary(j+1)/(xpt-x_cheby(j+1));
            numerator +=  t*functionvals[j];
            denominator +=  t;
         }
     }
    if (!var1){
         *output = numerator/denominator;
    }

}












void Mesh::InitDomain(const int Testcase, int *fixedDomain,int *fixedDisc, int *NelemX,int *NelemY,bool *PeriodicBD_X, bool*PeriodicBD_Y, dfloat *xL, dfloat *xR, dfloat *yL, dfloat *yR){

 switch(Testcase){
case 1: {    // convergence test, periodic
    *fixedDomain = 1 ;
    *fixedDisc = 0 ;
    *xL=0.0;
    *xR=1.0;
    *yL=0.0;
    *yR=1.0;
    break;}

case 2: {    // WELL BALANCED (CARTESIAN 20x20)
    *fixedDomain = 0 ;
    *fixedDisc = 1 ;
    *NelemX=20;
    *NelemY=20;
    *PeriodicBD_X=0;
    *PeriodicBD_Y=0;
    break;}
case 3:{    // Steeper Dam Break To Test Shock Capturing
    *fixedDomain = 1 ;
    *fixedDisc = 0;

    *xL=-5.0;
    *xR=5.0;
    *yL=-5.0;
    *yR=5.0;
    break;}
case 4: {    // WELL BALANCED (CARTESIAN 20x20)
    *fixedDomain = 0 ;
    *fixedDisc = 1 ;
    *NelemX=4;
    *NelemY=4;
    *PeriodicBD_X=0;
    *PeriodicBD_Y=0;
    break;}
case 30:{    // Steeper Dam Break To Test Shock Capturing
    *fixedDomain = 1 ;
    *fixedDisc = 0;

    *xL=-5.0;
    *xR=5.0;
    *yL=-5.0;
    *yR=5.0;
    break;}
case 31:     // Two-dimensional oscillating lake  (Xing_PosPres paper, 6.8)
    {
    *fixedDomain = 1 ;
    *fixedDisc = 0;

        *xL=-2.0;
        *xR=2.0;
        *yL=-2.0;
        *yR=2.0;


    break;

    }
case 32:     // Three Mound (4.6)
    {
    *fixedDomain = 1 ;
    *fixedDisc = 1;
    *xL=0.0;
    *xR=75.0;
    *yL=0.0;
    *yR=30.0;
    *NelemX=150;
    *NelemY=100;
    *PeriodicBD_X=0;
    *PeriodicBD_Y=0;


    break;

    }

case 33:     // Dam Break Three Mound (4.6)
    {

    *fixedDomain = 1 ;
    *fixedDisc = 1;
    *xL=0.0;
    *xR=75.0;
    *yL=0.0;
    *yR=30.0;
    *NelemX=150;
    *NelemY=100;
    *PeriodicBD_X=0;
    *PeriodicBD_Y=0;


    break;

    }


case 34:{     // 2D Solitary Wave Runup and Run-down
    *fixedDomain = 1 ;
    *fixedDisc = 0;

    *xL=0.0;
    *xR=25.0;
    *yL=-15.0;
    *yR=15.0;


    break;}
case 35:{     // 1D Bowl
    *fixedDomain = 1 ;
    *fixedDisc = 0;

    *xL=-3000.0;
    *xR=3000.0;
    *yL=-3000.0;
    *yR=3000.0;


    break;}

default:{
    *fixedDomain = 0 ;
    *fixedDisc = 0 ;
    break;}
    }


  }



