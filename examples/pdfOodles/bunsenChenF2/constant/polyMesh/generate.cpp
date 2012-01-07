#include "infra.H"

#define NR 5

int main(int argc, char* argv[])
{

    double r1=6, r2=7.5, /*r3=145.5,*/r3=34.0, r4=145.5, L=6.0;
  int nL=15, nAx=180;
  double axGrad=10.;
  int nr[NR-1]={8,2,25,25};
  double Lin=100;
  int nin=50;

  Matrix I(3,3,0.0);
  I(0,0)=1.;
  I(1,1)=1.;
  I(2,2)=1.;
  
  Matrix dm[4]={
    I,
    turn(-90.),
    turn(-180.),
    turn(-270.)
  };

  Vector dL(400., 0, 0), dI(-Lin,0,0);

  Vector front[4*NR];

  for (int i=0;i<4;i++)
    {
      front[NR*i+0]=dm[i]*Vector(0.,L/2.,0.);
      front[NR*i+1]=dm[i]*Vector(0.,r1,0.);
      front[NR*i+2]=dm[i]*Vector(0.,r2,0.);
      front[NR*i+3]=dm[i]*Vector(0.,r3,0.);
      front[NR*i+4]=dm[i]*Vector(0.,r4,0.);
    }
  //Kern
  blocks.push_back
    (
     Block
     (
      front[0], front[0]+dL, 
      front[3*NR]+dL, front[3*NR],
      front[1*NR], front[1*NR]+dL, 
      front[2*NR]+dL, front[2*NR],
      nAx, nL, nL,
      axGrad
      )
     );
  patches["axial_outlet"].addFace(front[0]+dL, front[1*NR]+dL, front[2*NR]+dL, front[3*NR]+dL);
/*
  //Kern verlängert
  blocks.push_back
    (
     Block
     (
      front[0]+dI, front[0], 
      front[3*NR], front[3*NR]+dI,
      front[1*NR]+dI, front[1*NR], 
      front[2*NR], front[2*NR]+dI,
      nin, nL, nL
      )
     );
  patches["inlet"].addFace(front[0]+dI, front[3*NR]+dI, front[2*NR]+dI, front[1*NR]+dI);
  */
  patches["inlet"].addFace(front[0], front[3*NR], front[2*NR], front[1*NR]);

  for (int i=0;i<4;i++)
    {
      int b=NR*i;
      int e;
      if (i<3) e=NR*(i+1);
      else e=0;
      //edges.push_back(new CircularEdge(front[b+1]+dI, front[e+1]+dI));
      for (int j=1;j<NR;j++)
	{
	  edges.push_back(new CircularEdge(front[b+j], front[e+j]));
	  edges.push_back(new CircularEdge(front[b+j]+dL, front[e+j]+dL));
	  patches["axial_outlet"].addFace(front[b+j-1]+dL, front[e+j-1]+dL, front[e+j]+dL, front[b+j]+dL);  
	}
      blocks.push_back
	(
	 Block
	 (
	  front[b+0], front[b+0]+dL, 
	  front[e+0]+dL, front[e+0],
	  front[b+1], front[b+1]+dL, 
	  front[e+1]+dL, front[e+1],
	  nAx, nL, nr[0],
	  axGrad
	  )
	 );
      /* * * Verlängerung * * */
/*
      blocks.push_back
	(
	 Block
	 (
	  front[b+0]+dI, front[b+0], 
	  front[e+0], front[e+0]+dI,
	  front[b+1]+dI, front[b+1], 
	  front[e+1], front[e+1]+dI,
	  nin, nL, nr[0]
	  )
	 );
      patches["walls"].addFace(front[b+1], front[b+1]+dI, front[e+1]+dI, front[e+1]);
      patches["inlet"].addFace(front[b+0]+dI, front[e+0]+dI, front[e+1]+dI, front[b+1]+dI);
*/
      patches["inlet"].addFace(front[b+0], front[e+0], front[e+1], front[b+1]);
      /* * * * * */
      blocks.push_back
	(
	 Block
	 (
	  front[b+1], front[b+1]+dL, 
	  front[e+1]+dL, front[e+1],
	  front[b+2], front[b+2]+dL, 
	  front[e+2]+dL, front[e+2],
	  nAx, nL, nr[1],
	  axGrad
	  )
	 );
      blocks.push_back
	(
	 Block
	 (
	  front[b+2], front[b+2]+dL, 
	  front[e+2]+dL, front[e+2],
	  front[b+3], front[b+3]+dL, 
	  front[e+3]+dL, front[e+3],
	  nAx, nL, nr[2],
	  axGrad, 1, 1.
	  )
	 );
      blocks.push_back
	(
	 Block
	 (
	  front[b+3], front[b+3]+dL, 
	  front[e+3]+dL, front[e+3],
	  front[b+4], front[b+4]+dL, 
	  front[e+4]+dL, front[e+4],
	  nAx, nL, nr[3],
	  axGrad, 1, 8.
	  )
	 );
      //patches["inlet"].addFace(front[b+0], front[e+0], front[e+1], front[b+1]);
      patches["walls"].addFace(front[b+1], front[e+1], front[e+2], front[b+2]);
      patches["coflow"].addFace(front[b+2], front[e+2], front[e+3], front[b+3]);
      patches["coflow_air"].addFace(front[b+3], front[e+3], front[e+4], front[b+4]);
      patches["tangential_outlet"].addFace(front[b+4], front[e+4], front[e+4]+dL, front[b+4]+dL);
    }
 
 
 numberVertices();

 std::ofstream f("blockMeshDict");
 f<<"FoamFile"
   <<std::endl<<"{"<<std::endl
   <<"version         2.0;"<<std::endl
   <<"format          ascii;"<<std::endl
   <<"root            \"\";"<<std::endl
   <<"case            \"\";"<<std::endl
   <<"instance        \"\";"<<std::endl
   <<"local           \"\";"<<std::endl
   <<"class           dictionary;"<<std::endl
   <<"object          blockMeshDict;"<<std::endl
   <<"}"<<std::endl
   <<"convertToMeters 0.001;"<<std::endl;

 f<<"vertices"<<std::endl<<"("<<std::endl;
 for 
  (
   std::map<Vector,int>::iterator it=vertexIDs.begin();
   it!=vertexIDs.end();
   it++
  )
  f<<" ("<<it->first<<")"<<std::endl;
  f<<");"<<std::endl;

 f<<"blocks"<<std::endl<<"("<<std::endl;
 for (std::vector<Block>::iterator it=blocks.begin();
       it!=blocks.end();it++)
 it->write(f);
 f<<");"<<std::endl;

 f<<"edges"<<std::endl<<"("<<std::endl;
 for (std::vector<ArcEdge*>::iterator it=edges.begin();
       it!=edges.end();it++)
 (*it)->write(f);
 f<<");"<<std::endl;

 f<<"patches"<<std::endl<<"("<<std::endl;
 for (std::map<std::string, Patch>::iterator it=patches.begin();
       it!=patches.end();it++)
   it->second.write(f, it->first);
 f<<");"<<std::endl;
 /*
 f<<"mergePatchPairs"<<std::endl<<"("<<std::endl;
 f<<" (wall_low_ax ax_out)"<<std::endl;
 f<<" (wall_low_tan tan_out)"<<std::endl;
 f<<");"<<std::endl;
 */
 return 0;
} 
