using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace OpenSees_NET_x64.EX1
{
    class Program
    {
        static void Main(string[] args)
        {

            var id = new OpenSees.IDWrapper(10);
            id[1] = 2;
            var mat = new OpenSees.MatrixWrapper(2, 2);
            mat[1, 1] = 3;
            var theDomain = new OpenSees.Components.DomainWrapper();
            var node1 = new OpenSees.Components.NodeWrapper(1, 2, 0, 0);
            var node2 = new OpenSees.Components.NodeWrapper(2, 2, 144.0, 0);
            var node3 = new OpenSees.Components.NodeWrapper(3, 2, 168.0, 0);
            var node4 = new OpenSees.Components.NodeWrapper(4, 2, 72.0, 96.0);

            theDomain.AddNode(new OpenSees.Components.NodeWrapper[] { node1, node2, node3, node4 });

            var theMaterial = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(1, 3000, 0);
            var truss1 = new OpenSees.Elements.TrussWrapper(1, 2, 1, 4, theMaterial, 10.0, 0, 0, 0);
            var truss2 = new OpenSees.Elements.TrussWrapper(2, 2, 2, 4, theMaterial, 5.0, 0, 0, 0);
            var truss3 = new OpenSees.Elements.TrussWrapper(3, 2, 3, 4, theMaterial, 5.0, 0, 0, 0);
            theDomain.AddElement(new OpenSees.Elements.ElementWrapper[] { truss1, truss2, truss3 });

            var sp1 = new OpenSees.Components.Constraints.SP_ConstraintWrapper(1, 0, 0.0, true);
            var sp2 = new OpenSees.Components.Constraints.SP_ConstraintWrapper(1, 1, 0.0, true);
            var sp3 = new OpenSees.Components.Constraints.SP_ConstraintWrapper(2, 0, 0.0, true);
            var sp4 = new OpenSees.Components.Constraints.SP_ConstraintWrapper(2, 1, 0.0, true);
            var sp5 = new OpenSees.Components.Constraints.SP_ConstraintWrapper(3, 0, 0.0, true);
            var sp6 = new OpenSees.Components.Constraints.SP_ConstraintWrapper(3, 1, 0.0, true);
            theDomain.AddSP_Constraint(new OpenSees.Components.Constraints.SP_ConstraintWrapper[] { sp1, sp2, sp3, sp4, sp5, sp6 });

            var theSeries = new OpenSees.Components.Timeseries.LinearSeriesWrapper();
            var theLoadPattern = new OpenSees.Components.LoadPatterns.LoadPatternWrapper(1);
            theLoadPattern.SetTimeSeries(theSeries);
            theDomain.AddLoadPattern(theLoadPattern);

            var theLoadValues = new OpenSees.VectorWrapper(2);
            theLoadValues[0] = 100;
            theLoadValues[1] = -50;

            var nodalLoad = new OpenSees.Components.Loads.NodalLoadWrapper(1, 4, theLoadValues, false);
            theDomain.AddNodalLoad(nodalLoad, 1);

            var theModel = new OpenSees.AnalysisModelWrapper();
            var theSolnAlgo = new OpenSees.Algorithms.LinearWrapper();
            var theIntegrator = new OpenSees.Integrators.Static.LoadControlWrapper(1.0, 1, 1.0, 1.0);
            var theHandler = new OpenSees.Handlers.PlainHandlerWrapper();
            var theRCM = new OpenSees.GraphNumberers.RCMWrapper(false);
            var theNumberer = new OpenSees.Numberers.DOF_NumbererWrapper(theRCM);
            var theSolver = new OpenSees.Systems.Linears.BandSPDLinLapackSolverWrapper();
            var theSOE = new OpenSees.Systems.Linears.BandSPDLinSOEWrapper(theSolver);
            var theAnalysis = new OpenSees.Analysis.StaticAnalysisWrapper(
                theDomain,
                theHandler,
                theNumberer,
                theModel,
                theSolnAlgo,
                theSOE,
                theIntegrator,
                null);

            int numSteps = 1;
            var ret = theAnalysis.Analyze(numSteps);

            theDomain.CalculateNodalReactions(0);
            var disp = node4.GetCommitDisp();
            Console.WriteLine($"node disps : {{{string.Join(" ,", disp)}}}");
            Console.ReadKey();
        }
    }
}
