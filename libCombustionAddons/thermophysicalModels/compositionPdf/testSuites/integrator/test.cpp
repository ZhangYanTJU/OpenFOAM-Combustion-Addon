#include "fvCFD.H"
#include "chemistryTable.H"
#include "integrator.H"
#include "beta.H"

#include <vector>
#include <fstream>

double fac=1.0;
beta globybetapdf("y");

sumOfIntegrateableFunctions<scalar> pdfs(double mean, double var)
{
    globybetapdf.setParameters(mean, var*mean*(1.0-mean)+mean*mean);
    return fac * globybetapdf;
}

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
            20,            
            20
        };
    double start[]=
        {
            0,
            0,
            0
        };
    double step[]=
        {
            1.0/double(res[0]-1), 
            1.0/double(res[1]-1),
            1.0/double(res[2]-1)
        };

    int dim=3;

    label n=0;
    HashTable<label,word> pvdesc;
    pvdesc.insert("x",n++);
    pvdesc.insert("y",n++);
    pvdesc.insert("z",n++);
    n=0;
    HashTable<label,word> cdesc;
    cdesc.insert("c",n++);
    
    chemistryTable tab
        (
            dim, res, start, step,
            cdesc, pvdesc
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
            for 
                (
                    i[2]=0;
                    i[2]<tab.nElements(2);
                    i[2]++
                )
            {
                double x=tab.valueAt(0, i[0]);
                double y=tab.valueAt(1, i[1]);
                double z=tab.valueAt(2, i[2]);
                chemicalSystemState V(cdesc);
                V[0]=
                    /*1.0 +*/ 5*x*x+7*y*y*y+9*z
                    ;
                tab.access(i)=V;

            }
    Info<<"Table set up: c=5*x*x+7*y*y*y+9*z"<<endl;




    Info<<endl<<endl<<"== Testing plain integration =="<<endl<<endl;

    integrationParameter x;
    integrationParameter constants;
    scalar Ic;

    x.insert("x", 0.0);

    constants.insert("y", 0.0);
    constants.insert("z", 0.0);
    Ic=integrator<scalar>::integrate
        (
            tab.entry("c"),
            x, constants
        );
    Info<<"(y=0, z=0) int from x=0 to 1: I="<<Ic<<endl;

    constants.clear();
    constants.insert("y", 0.5);
    constants.insert("z", 0.5);
    Ic=integrator<scalar>::integrate
        (
            tab.entry("c"),
            x, constants
        );
    Info<<"(y=0.5, z=0.5) int from x=0 to 1: I="<<Ic<<endl;

    
    x.insert("y", 0.0);
    x.insert("z", 0.0);
    constants.clear();
    Ic=integrator<scalar>::integrate
        (
            tab.entry("c"),
            x, constants
        );
    Info<<"int from x=0 to 1 (int from y=0 to 1 (int from z=0 to 1)): I="<<Ic<<endl;




    Info<<endl<<endl<<"== Testing combined integration =="<<endl<<endl;

    x.clear();
    constants.clear();

    x.insert("x", 0.0);
    constants.insert("y", 0.0);
    constants.insert("z", 0.0);

    Ic=integrator<scalar>::integrate
        (
            singleParameter("x")*tab.entry("c"),
            x, constants
        );
    Info<<"(y=0, z=0) int from x=0 to 1 (x*F): I="<<Ic<<endl;
    
    constants.clear();
    constants.insert("y", 0.5);
    constants.insert("z", 0.5);
    Ic=integrator<scalar>::integrate
        (
            singleParameter("x")*singleParameter("x")*tab.entry("c"),
            x, constants
        );
    Info<<"(y=0.5, z=0.5) int from x=0 to 1 (x*x*F): I="<<Ic<<endl;

    
    constants.clear();
    x.insert("y", 0.0);
    x.insert("z", 0.0);
    Ic=integrator<scalar>::integrate
        (
            singleParameter("x")*singleParameter("x")*tab.entry("c"),
            x, constants
        );
    Info<<"int from x=0 to 1 (int from y=0 to 1 (int from z=0 to 1)) (x*x*F): I="<<Ic<<endl;

    Ic=integrator<scalar>::integrate
        (
            singleParameter("x")*singleParameter("x")*tab.entry("c")*tab.entry("c"),
            x, constants
        );
    Info<<"int from x=0 to 1 (int from y=0 to 1 (int from z=0 to 1)) (x*x*F*F): I="<<Ic<<endl;



    Info<<endl<<endl<<"== Testing smooth PDF integration =="<<endl<<endl;

    beta smoothbetapdf("x");
    smoothbetapdf.setParameters(0.285714, 0.107143);
    beta smoothbetapdf2("y");
    smoothbetapdf2.setParameters(0.285714, 0.107143);
    beta smoothbetapdf3("z");
    smoothbetapdf3.setParameters(0.5, 0.3);

    x.clear();
    constants.clear();

    x.insert("x", 0.0);
    constants.insert("y", 0.0);
    constants.insert("z", 0.0);

    Ic=integrator<scalar>::integrate
        (
            tab.entry("c") * smoothbetapdf,
            x, constants
        );
    Info<<"(y=0, z=0) int from x=0 to 1 (F*betapdf): I="<<Ic<<endl;
 
    x.clear();
    constants.clear();
    x.insert("x", 0.0);
    x.insert("y", 0.0);
    constants.insert("z", 0.0);

    Ic=integrator<scalar>::integrate
        (
            tab.entry("c") * smoothbetapdf * smoothbetapdf2,
            x, constants
        );
    Info<<"(y=0, z=0) int from x=0 to 1 (F*betapdf*betapdf2): I="<<Ic<<endl;

    Ic=integrator<scalar>::integrate
        (
            singleParameter("x") * singleParameter("y")
              * tab.entry("c") * smoothbetapdf * smoothbetapdf2,
            x, constants
        );
    Info<<"(y=0, z=0) int from x=0 to 1 (x*y*F*betapdf*betapdf2): I="<<Ic<<endl;

    x.clear();
    constants.clear();
    x.insert("x", 0.0);
    x.insert("y", 0.0);
    x.insert("z", 0.0);
    Ic=integrator<scalar>::integrate
        (
            singleParameter("x") * singleParameter("y")
              * tab.entry("c") * smoothbetapdf * smoothbetapdf2 * smoothbetapdf3,
            x, constants
        );
    Info<<"int from x=0 to 1 (int from y=0 to 1 (int from z=0 to 1 (x*y*F*betapdf*betapdf2*betapdf3))): I="
        <<Ic<<endl;

    Ic=integrator<scalar>::integrate
        (
            smoothbetapdf * smoothbetapdf2 * smoothbetapdf3,
            x, constants
        );
    Info<<"int from x=0 to 1 (int from y=0 to 1 (int from z=0 to 1 (betapdf*betapdf2*betapdf3))): I="
        <<Ic<<endl;

    Info<<endl<<endl<<"== Testing PDF integration with boundary singularities =="<<endl<<endl;

    beta betapdf("x");
    betapdf.setParameters(0.5, 0.375);

    x.clear();
    constants.clear();

    x.insert("x", 0.0);
    constants.insert("y", 0.0);
    constants.insert("z", 0.0);

    Ic=integrator<scalar>::integrate
        (
            tab.entry("c") * betapdf,
            x, constants
        );
    Info<<"(y=0, z=0) int from x=0 to 1 (F*betapdf): I="<<Ic<<endl;

    Info<<endl<<endl<<"== Testing PDF integration with boundary dirac peaks =="<<endl<<endl;

    betapdf.setParameters(0.5, 0.49995);

    x.clear();
    constants.clear();

    x.insert("x", 0.0);
    constants.insert("y", 0.0);
    constants.insert("z", 0.0);

    Ic=integrator<scalar>::integrate
        (
            tab.entry("c") * betapdf,
            x, constants
        );
    Info<<"(y=0, z=0) int from x=0 to 1 (F*betapdf): I="<<Ic<<endl;

    Info<<endl<<endl<<"== Testing PDF integration with interior dirac peak =="<<endl<<endl;

    beta ybetapdf("y");
    ybetapdf.setParameters(0.5, 0.250125);

    x.clear();
    constants.clear();

    x.insert("y", 0.0);
    constants.insert("x", 0.0);
    constants.insert("z", 0.0);

    Ic=integrator<scalar>::integrate
        (
            tab.entry("c") * ybetapdf,
            x, constants
        );
    Info<<"(x=0, z=0) int from y=0 to 1 (F*ybetapdf): I="<<Ic<<endl;

    Info<<endl<<endl<<"== Testing combined PDF integration smooth and with dirac peaks =="<<endl<<endl;

    x.clear();
    constants.clear();

    x.insert("x", 0.0);
    x.insert("y", 0.0);
    x.insert("z", 0.0);

    Ic=integrator<scalar>::integrate
        (
            tab.entry("c") * betapdf * ybetapdf * smoothbetapdf3,
            x, constants
        );
    Info<<"int from x=0 to 1 (int from y=0 to 1 ( int from z=0 to 1 (F*betapdf*ybetapdf*smoothbetapdf3))): I="
        <<Ic<<endl;




    Info<<endl<<endl<<"== Testing PDF integration loop =="<<endl<<endl;

double mean=0.5;
    for (int i=0; i<40; i++)
{
 double var=i*(1.0/(40.0-1.0));

    sumOfIntegrateableFunctions<scalar> copy=pdfs(mean, var);

    x.clear();
    constants.clear();

    x.insert("y", 0.0);
    constants.insert("x", 0.0);
    constants.insert("z", 0.0);

    Ic=integrator<scalar>::integrate
        (
            tab.entry("c") * /*ybetapdf*/ copy,
            x, constants
        );
    Info<<"I("<<mean<<","<<var<<")="<<Ic<<endl;
}




    return 0;
}
