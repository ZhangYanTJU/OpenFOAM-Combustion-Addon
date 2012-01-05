#include "fvCFD.H"
#include "MultidimensionalLookupTable.H"

#include <vector>
#include <fstream>


class Gnuplot
{
public:
  FILE *command;

  Gnuplot(string title)
  {
      command = popen("gnuplot -persist","w");
      //command=stdout;
    fprintf(command, "set title '%s'\n", title.c_str());
  }

  ~Gnuplot()
  {
    pclose(command);
  }

    void executeCommand(const string&cmd)
        {
            fprintf(command, "%s\n", cmd.c_str());
        }

    void startPlot(const std::vector<string>& curvetitles, string addplots="")
        {
            fprintf(command, "plot ");
            if (addplots!="")
                fprintf(command, "%s, ", addplots.c_str());
            int i=0;
            for (std::vector<string>::const_iterator it=curvetitles.begin();
                 it!=curvetitles.end(); it++)
            {
                if (it==curvetitles.begin())
                    fprintf(command, "'-' ");
                else
                    fprintf(command, "'' ");

                fprintf(command, "w p t '%s'", it->c_str());

                if (it!=curvetitles.end()-1)
                    fprintf(command, ",");
                else
                    fprintf(command, "\n");
                i++;
            }
        }


    void addDataPoint(double x, double y)
        {
            fprintf(command, "%lf %lf\n", x, y);
        }

    void finishPlot()
        {
            fprintf(command, "e\n");
        }
    void breakDataset()
        {
            fprintf(command, "\n\n");
        }
};



int main(int argc, char*argv[])
{
    // create tab
    int res[]=
        {
            30,
            20
#ifdef dim3D
            ,
            20
#endif
        };
    double start[]=
        {
            -0.25,
            -0.5
#ifdef dim3D
            ,0
#endif
        };
    double step[]=
        {
            2.0/double(res[0]-1), 
            2.0/double(res[1]-1)
#ifdef dim3D
            , 1.0/double(res[2]-1)
#endif
        };

    int dim=
#ifdef dim3D
        3;
#else
        2;
#endif
  
    MultidimensionalLookupTable<scalar> tab
        (
            dim, res, start, step
        );

    // fill
    int i[dim];

    for 
        (
            i[0]=0;
            i[0]<tab.nElements(0);
            i[0]++
        )
        for 
            (
                i[1]=0;
                i[1]<tab.nElements(1);
                i[1]++
            )
#ifdef dim3D
            for 
                (
                    i[2]=0;
                    i[2]<tab.nElements(2);
                    i[2]++
                )
#endif
            {
                double x=tab.valueAt(0, i[0]);
                double y=tab.valueAt(1, i[1]);
#ifdef dim3D
                double z=tab.valueAt(2, i[2]);
#endif
                tab.access(i)=
                    5*x+7*y
#ifdef dim3D
                    +9*z
#endif
                    ;
            }

    Gnuplot gplt("testplot");
    gplt.executeCommand("set key top left");

    std::vector<string> titles;
    titles.push_back("slice diagonal");
    titles.push_back("slice x");
    titles.push_back("slice y");
#ifdef dim3D
    titles.push_back("slice z");
#endif

#ifdef dim3D
    gplt.startPlot(titles, "5*x, 7*x, 9*x, (5+7+9)*x");
#else
    gplt.startPlot(titles, "5*x, 7*x, (5+7)*x");
#endif

    std::ofstream f("testplot");
    label np=200;
    for (int i=0;i<np;i++)
    {
        double x=-0.5+(2.0/double(np-1))*double(i);
        scalarField xv(dim, x);
        double y=tab.lookup(xv);
        gplt.addDataPoint(x, y);
        f<<x<<" "<<x<<" "<<y<<std::endl;
    }
    gplt.finishPlot();

    for (int i=0;i<np;i++)
    {
        double x=-0.5+(2.0/double(np-1))*double(i);
        scalarField xv(dim);
        xv[0]=x;
        xv[1]=0.0;
#ifdef dim3D
        xv[2]=0.0;
#endif
        scalar y=tab.lookup(xv);
        gplt.addDataPoint(x, y);
    }
    gplt.finishPlot();

    for (int i=0;i<np;i++)
    {
        double x=-0.5+(2.0/double(np-1))*double(i);
        scalarField xv(dim);
        xv[0]=0.0;
        xv[1]=x;
#ifdef dim3D
        xv[2]=0.0;
#endif
        scalar y=tab.lookup(xv);
        gplt.addDataPoint(x, y);
    }
    gplt.finishPlot();

#ifdef dim3D
    for (int i=0;i<np;i++)
    {
        double x=-0.5+(2.0/double(np-1))*double(i);
        scalarField xv(dim);
        xv[0]=0.0;
        xv[1]=0.0;
        xv[2]=x;
        scalar y=tab.lookup(xv);
        gplt.addDataPoint(x, y);
    }
    gplt.finishPlot();
#endif

/*
    List<scalar> c(2, 1.0);
    Info<<"Reverse lookup:"<<endl;
    Info<<"x(F=12,y=1) = "<<tab.reverseLookup(0, c, 0, 12)<<endl;
*/

    return 0;
}
