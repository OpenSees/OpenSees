using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TestTclInterp
{
    class Program
    {
        static void Main(string[] args)
        {
            var tclInterp = new OpenSees.Tcl.TclWrapper();
            tclInterp.Init();
            var domain = tclInterp.GetActiveDomain();
            tclInterp.TclEvalFile("model.tcl");
            tclInterp.TclEval("print node 1");
            var node1 = domain.GetNode(tag:1);
            var disp = node1.GetCommitDisp();
            Console.WriteLine(string.Join(",", disp));
            Console.ReadKey();
        }
    }
}
