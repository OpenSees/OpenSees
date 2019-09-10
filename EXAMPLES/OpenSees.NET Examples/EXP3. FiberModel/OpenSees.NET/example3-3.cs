/// <summary>
        /// Example 3.2 : Transeint Analysis analysis
        /// </summary>
        public static void Exp3_3()
        {
            Console.WriteLine("Example 3.3 : Transeint Analysis analysis");
            var theDomain = new OpenSees.Components.DomainWrapper();

            double width = 360;
            double height = 144;
            var node1 = new OpenSees.Components.NodeWrapper(1, 3, 0, 0);
            var node2 = new OpenSees.Components.NodeWrapper(2, 3, width, 0);
            var node3 = new OpenSees.Components.NodeWrapper(3, 3, 0, height);
            var node4 = new OpenSees.Components.NodeWrapper(4, 3, width, height);
            var node5 = new OpenSees.Components.NodeWrapper(5, 3, 0, 2 * height);
            var node6 = new OpenSees.Components.NodeWrapper(6, 3, width, 2 * height);
            var node7 = new OpenSees.Components.NodeWrapper(7, 3, 0, 3 * height);
            var node8 = new OpenSees.Components.NodeWrapper(8, 3, width, 3 * height);
            theDomain.AddNode(new OpenSees.Components.NodeWrapper[] { node1, node2, node3,
                    node4, node5, node6, node7, node8 });

            var sp1 = new OpenSees.Components.Constraints.SP_ConstraintWrapper(1, 0, 0.0, true);
            var sp2 = new OpenSees.Components.Constraints.SP_ConstraintWrapper(1, 1, 0.0, true);
            var sp3 = new OpenSees.Components.Constraints.SP_ConstraintWrapper(1, 2, 0.0, true);
            var sp4 = new OpenSees.Components.Constraints.SP_ConstraintWrapper(2, 0, 0.0, true);
            var sp5 = new OpenSees.Components.Constraints.SP_ConstraintWrapper(2, 1, 0.0, true);
            var sp6 = new OpenSees.Components.Constraints.SP_ConstraintWrapper(2, 2, 0.0, true);
            theDomain.AddSP_Constraint(new OpenSees.Components.Constraints.SP_ConstraintWrapper[] { sp1, sp2, sp3,
                    sp4, sp5, sp6 });

            var corcrt = new OpenSees.Materials.Uniaxials.Concrete01Wrapper(1, -6.0, -0.004, -5.0, -0.014);
            var covcrt = new OpenSees.Materials.Uniaxials.Concrete01Wrapper(2, -5.0, -0.002, 0.0, -0.006);

            double fy = 60.0;
            double E = 30000.0;
            var stl = new OpenSees.Materials.Uniaxials.Steel01Wrapper(3, fy, E, 0.01);

            double colWidth = 15.0;
            double colDepth = 24.0;
            double cover = 1.5;
            double As = 0.6;

            var y1 = colDepth / 2;
            var z1 = colWidth / 2;

            var repres = new OpenSees.Materials.Sections.Repres.FiberSectionReprWrapper(1,
                new OpenSees.Materials.Sections.Repres.PatchWrapper[]
                {
                        new OpenSees.Materials.Sections.Repres.QuadPatchWrapper(1,10,1,new OpenSees.MatrixWrapper(new double[,]{{ (y1 - cover), (z1 - cover) },{ -(y1 - cover), (z1 - cover) },{ -(y1 - cover), -(z1 - cover) },{ (y1 - cover), -(z1 - cover) } })),
                        new OpenSees.Materials.Sections.Repres.QuadPatchWrapper(2,10,1,new OpenSees.MatrixWrapper(new double[,]{{ (y1), (z1) },{ -(y1), (z1) },{ -(y1), (z1 - cover) },{ (y1), (z1 - cover) } })),
                        new OpenSees.Materials.Sections.Repres.QuadPatchWrapper(2,10,1,new OpenSees.MatrixWrapper(new double[,]{{ (y1), -(z1 - cover) },{ -(y1), -(z1 - cover) },{ -(y1), -(z1) },{(y1), -(z1) } })),
                        new OpenSees.Materials.Sections.Repres.QuadPatchWrapper(2,2,1,new OpenSees.MatrixWrapper(new double[,]{{ (y1), (z1-cover) },{ (y1 - cover), (z1 - cover) },{ (y1 - cover), -(z1 - cover) },{ (y1), -(z1-cover) } })),
                        new OpenSees.Materials.Sections.Repres.QuadPatchWrapper(2,2,1,new OpenSees.MatrixWrapper(new double[,]{{ -(y1-cover), (z1-cover) },{ -(y1), (z1 - cover) },{ -(y1), -(z1 - cover) },{ -(y1-cover), -(z1-cover) } })),
                },
                new OpenSees.Materials.Sections.Repres.ReinfLayerWrapper[]
                {
                        new OpenSees.Materials.Sections.Repres.StraightReinfLayerWrapper(3,3,As,new OpenSees.VectorWrapper(new double[]{ -y1+cover,z1-cover }), new OpenSees.VectorWrapper(new double[]{ -y1+cover,-z1+cover})),
                        new OpenSees.Materials.Sections.Repres.StraightReinfLayerWrapper(3,2,As,new OpenSees.VectorWrapper(new double[]{ 0,z1-cover }), new OpenSees.VectorWrapper(new double[]{ 0,-z1+cover})),
                        new OpenSees.Materials.Sections.Repres.StraightReinfLayerWrapper(3,3,As,new OpenSees.VectorWrapper(new double[]{ y1-cover,z1-cover }), new OpenSees.VectorWrapper(new double[]{ y1-cover,-z1+cover})),
                });

            var section = new OpenSees.Materials.Sections.FiberSection2dWrapper(1, repres,
                new Dictionary<int, OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper>() {
                        {1,corcrt },
                        {2,covcrt },
                        {3,stl }
                });

            var goeTransCol = new OpenSees.Elements.CrdTransfs.CorotCrdTransf2dWrapper();
            var np = 5;
            var col1 = new OpenSees.Elements.ForceBeamColumn2dWrapper(1, 1, 3, np, section, new OpenSees.Elements.BeamIntegrations.BeamIntegrationWrapper(OpenSees.Elements.BeamIntegrations.BeamIntegrationType.Lobatto), goeTransCol);
            var col2 = new OpenSees.Elements.ForceBeamColumn2dWrapper(2, 2, 4, np, section, new OpenSees.Elements.BeamIntegrations.BeamIntegrationWrapper(OpenSees.Elements.BeamIntegrations.BeamIntegrationType.Lobatto), goeTransCol);
            var col3 = new OpenSees.Elements.ForceBeamColumn2dWrapper(5, 3, 5, np, section, new OpenSees.Elements.BeamIntegrations.BeamIntegrationWrapper(OpenSees.Elements.BeamIntegrations.BeamIntegrationType.Lobatto), goeTransCol);
            var col4 = new OpenSees.Elements.ForceBeamColumn2dWrapper(6, 4, 6, np, section, new OpenSees.Elements.BeamIntegrations.BeamIntegrationWrapper(OpenSees.Elements.BeamIntegrations.BeamIntegrationType.Lobatto), goeTransCol);
            var col5 = new OpenSees.Elements.ForceBeamColumn2dWrapper(8, 5, 7, np, section, new OpenSees.Elements.BeamIntegrations.BeamIntegrationWrapper(OpenSees.Elements.BeamIntegrations.BeamIntegrationType.Lobatto), goeTransCol);
            var col6 = new OpenSees.Elements.ForceBeamColumn2dWrapper(9, 6, 8, np, section, new OpenSees.Elements.BeamIntegrations.BeamIntegrationWrapper(OpenSees.Elements.BeamIntegrations.BeamIntegrationType.Lobatto), goeTransCol);

            var goeTransBeam = new OpenSees.Elements.CrdTransfs.LinearCrdTransf2dWrapper();
            var beam1 = new OpenSees.Elements.ElasticBeam2dWrapper(3, 360, 4030, 8640, 3, 4, goeTransBeam, 0, 0, 0, 0);
            var beam2 = new OpenSees.Elements.ElasticBeam2dWrapper(7, 360, 4030, 8640, 5, 6, goeTransBeam, 0, 0, 0, 0);
            var beam3 = new OpenSees.Elements.ElasticBeam2dWrapper(10, 360, 4030, 8640, 7, 8, goeTransBeam, 0, 0, 0, 0);

            theDomain.AddElement(new OpenSees.Elements.ElementWrapper[] { col1, col2, col3, col4, col5, col6, beam1, beam2, beam3 });

            var theSeries = new OpenSees.Components.Timeseries.LinearSeriesWrapper();
            var gravityLoadPattern = new OpenSees.Components.LoadPatterns.LoadPatternWrapper(1);
            gravityLoadPattern.SetTimeSeries(theSeries);
            theDomain.AddLoadPattern(gravityLoadPattern);

            var gravityLoadValues = new OpenSees.VectorWrapper(new double[] { 0, -180, 0 });

            var nl1 = new OpenSees.Components.Loads.NodalLoadWrapper(1, 3, gravityLoadValues, false);
            var nl2 = new OpenSees.Components.Loads.NodalLoadWrapper(2, 4, gravityLoadValues, false);
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper[] { nl1, nl2 }, 1);
            {
                var theModel = new OpenSees.AnalysisModelWrapper();
                var theSolnAlgo = new OpenSees.Algorithms.NewtonRaphsonWrapper();
                var theIntegrator = new OpenSees.Integrators.Static.LoadControlWrapper(0.1, 1, 0.1, 0.1);
                var theHandler = new OpenSees.Handlers.TransformationConstraintHandlerWrapper();
                var theRCM = new OpenSees.GraphNumberers.RCMWrapper(false);
                var theNumberer = new OpenSees.Numberers.DOF_NumbererWrapper(theRCM);
                var theSolver = new OpenSees.Systems.Linears.ProfileSPDLinDirectSolverWrapper();
                var theSOE = new OpenSees.Systems.Linears.ProfileSPDLinSOEWrapper(theSolver);
                var theTest = new OpenSees.ConvergenceTests.CTestNormDispIncrWrapper(1e-12, 10, 3, 2, 1.0e10);
                var theAnalysis = new OpenSees.Analysis.StaticAnalysisWrapper(
                    theDomain,
                    theHandler,
                    theNumberer,
                    theModel,
                    theSolnAlgo,
                    theSOE,
                    theIntegrator,
                    theTest);

                theAnalysis.Analyze(10);
            }

            node3.PrintSelf(0);
            node4.PrintSelf(0);

            Console.WriteLine("Gravity load analysis completed");


            Console.WriteLine("Transeint Analysis analysis");
            {
                theDomain.SetLoadConst();
                theDomain.SetCurrentTime(0);
                theDomain.SetCommittedTime(0);
                theDomain.RemoveRecorders();

                double g = 386.4;
                double p = 180;
                double m = p / g;
                var mass = new double[,] { { m, 0, 0 }, { 0, m, 0 }, { 0, 0, 0 } };
                var massMat = new OpenSees.MatrixWrapper(mass);
                theDomain.SetMass(3, massMat);
                theDomain.SetMass(4, massMat);

                double dt = 0.02;
                var accelVector = GetFileDataSeries($@"{Environment.CurrentDirectory}\ARL360.at2", 8);
                var accelSeries = new OpenSees.Components.Timeseries.PathSeriesWrapper(accelVector, dt, g, false, false, 0);
                var groundMotion = new OpenSees.Components.GroundMotions.GroundMotionWrapper(null, null, accelSeries, null, dt, 1.0);
                var uniformExcitationPattern = new OpenSees.Components.LoadPatterns.UniformExcitationWrapper(2, groundMotion, 0, 0, 1);
                theDomain.AddLoadPattern(uniformExcitationPattern);

                theDomain.SetRayleighDampingFactors(0, 0, 0, 0.000625);

                {
                    var theModel = new OpenSees.AnalysisModelWrapper();
                    var theSolnAlgo = new OpenSees.Algorithms.NewtonRaphsonWrapper();
                    var theIntegrator = new OpenSees.Integrators.Transient.NewmarkWrapper(0.5, 0.25);
                    var theHandler = new OpenSees.Handlers.PlainHandlerWrapper();
                    var theRCM = new OpenSees.GraphNumberers.RCMWrapper(false);
                    var theNumberer = new OpenSees.Numberers.DOF_NumbererWrapper(theRCM);
                    var theSolver = new OpenSees.Systems.Linears.ProfileSPDLinDirectSolverWrapper();
                    var theSOE = new OpenSees.Systems.Linears.ProfileSPDLinSOEWrapper(theSolver);
                    var theTest = new OpenSees.ConvergenceTests.CTestNormDispIncrWrapper(1.0e-12, 10, 3, 2, 1e10);
                    var theAnalysis = new OpenSees.Analysis.DirectIntegrationAnalysisWrapper(
                        theDomain,
                        theHandler,
                        theNumberer,
                        theModel,
                        theSolnAlgo,
                        theSOE,
                        theIntegrator,
                        theTest);

                    var theEigenSOE = new OpenSees.Systems.Eigens.FullGenEigenSOEWrapper(new OpenSees.Systems.Eigens.FullGenEigenSolverWrapper(), theModel);
                    theAnalysis.SetEigenSOE(theEigenSOE);

                    // eigen analysis
                    {
                        // puts "eigen values at start of transient: [eigen 2]"
                        var ret = theAnalysis.Eigen(2, true, true);
                        if (ret == 0)
                        {
                            var eigens = theDomain.GetEigenValues();
                            Console.WriteLine($"eigen values at start of transient: {string.Join(" ", eigens)}");
                        }
                    }

                    // recorders
                    {
                        var opsstream = new OpenSees.Handlers.DataFileStreamWrapper(Environment.CurrentDirectory + $@"\Exp3_1_nodes_disp.txt");
                        var recorder = new OpenSees.Recorders.NodeRecorderWrapper(new OpenSees.IDWrapper(new int[]{ 0,1,2}), new OpenSees.IDWrapper(new int[] { 3,4 }), 0, "disp", theDomain, opsstream);
                        theDomain.AddRecorder(recorder);
                    }

                    {
                        var adt = 0.01;
                        var tfinal = 2000 * 0.01;
                        var tcurrent = theDomain.GetTime();
                        var ok = 0;
                        do
                        {
                            ok = theAnalysis.Analyze(1, adt);
                            if (ok != 0) // analysis failed
                            {
                                Console.WriteLine($"Transeint Analysis Failed!");
                                break;
                            }
                            tcurrent = theDomain.GetTime();
                        } while (tcurrent < tfinal);
                        if (ok == 0)
                            Console.WriteLine($"Transeint Analysis Successfull Done!");
                    }
                }

                node3.PrintSelf(0);
                Console.WriteLine("Transeint Analysis completed");
            }
        }


GetFileDataSeries(string filename
            , int skipline = 0, double factor = 1)
        {
            var sr = new System.IO.StreamReader(filename);
            for (var i = 0; i < skipline; i++)
                sr.ReadLine();
            var data = new List<double>();
            while (!sr.EndOfStream)
            {
                var lineData = sr.ReadLine().Split(new char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries).Select(i => double.Parse(i) * factor);
                if (lineData.Count() > 0)
                    data.AddRange(lineData);
            }
            return new OpenSees.VectorWrapper(data.ToArray());
        }