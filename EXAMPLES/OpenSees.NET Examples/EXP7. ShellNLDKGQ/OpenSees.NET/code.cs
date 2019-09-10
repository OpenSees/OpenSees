//http://www.luxinzheng.net/download/OpenSEES/Examples_of_NLDKGQ_element.htm
        public static void Exp7_ShellNLDKGQ()
        {
            var theDomain = new OpenSees.Components.DomainWrapper();
            var node1 = new OpenSees.Components.NodeWrapper(1, 6, 0.0, 1.0, 0.0);//node 1          0.   1.0   0.
            theDomain.AddNode(node1);

            var node2 = new OpenSees.Components.NodeWrapper(2, 6, 0.0, 0.5, 0.0);//node 2          0.   0.5   0.
            theDomain.AddNode(node2);

            var node3 = new OpenSees.Components.NodeWrapper(3, 6, 0.0, 0.0, 0.0);//node 3          0.    0.   0.
            theDomain.AddNode(node3);

            var node4 = new OpenSees.Components.NodeWrapper(4, 6, 1.0, 1.0, 0.0);//node 4         1.0    1.   0.
            theDomain.AddNode(node4);

            var node5 = new OpenSees.Components.NodeWrapper(5, 6, 1.0, 0.5, 0.0);//node 5         1.0   0.5   0.
            theDomain.AddNode(node5);

            var node6 = new OpenSees.Components.NodeWrapper(6, 6, 1.0, 0.0, 0.0);//node 6         1.0    0.   0.
            theDomain.AddNode(node6);

            var node7 = new OpenSees.Components.NodeWrapper(7, 6, 2.0, 1.0, 0.0);//node 7         2.0    1.   0.
            theDomain.AddNode(node7);

            var node8 = new OpenSees.Components.NodeWrapper(8, 6, 2.0, 0.5, 0.0);//node 8         2.0   0.5   0.
            theDomain.AddNode(node8);

            var node9 = new OpenSees.Components.NodeWrapper(9, 6, 2.0, 0.0, 0.0);//node 9         2.0    0.   0.
            theDomain.AddNode(node9);

            var node10 = new OpenSees.Components.NodeWrapper(10, 6, 3.0, 1.0, 0.0);//node 10        3.0    1.   0.
            theDomain.AddNode(node10);

            var node11 = new OpenSees.Components.NodeWrapper(11, 6, 3.0, 0.5, 0.0);//node 11        3.0   0.5   0.
            theDomain.AddNode(node11);

            var node12 = new OpenSees.Components.NodeWrapper(12, 6, 3.0, 0.0, 0.0);//node 12        3.0    0.   0.
            theDomain.AddNode(node12);

            var node13 = new OpenSees.Components.NodeWrapper(13, 6, 4.0, 1.0, 0.0);//node 13        4.0    1.   0.
            theDomain.AddNode(node13);

            var node14 = new OpenSees.Components.NodeWrapper(14, 6, 4.0, 0.5, 0.0);//node 14        4.0   0.5   0.
            theDomain.AddNode(node14);

            var node15 = new OpenSees.Components.NodeWrapper(15, 6, 4.0, 0.0, 0.0);//node 15        4.0    0.   0.
            theDomain.AddNode(node15);

            var node16 = new OpenSees.Components.NodeWrapper(16, 6, 5.0, 1.0, 0.0);//node 16         5.    1.   0.
            theDomain.AddNode(node16);

            var node17 = new OpenSees.Components.NodeWrapper(17, 6, 5.0, 0.5, 0.0);//node 17         5.   0.5   0.
            theDomain.AddNode(node17);

            var node18 = new OpenSees.Components.NodeWrapper(18, 6, 5.0, 0.0, 0.0);//node 18         5.    0.   0.
            theDomain.AddNode(node18);

            var node19 = new OpenSees.Components.NodeWrapper(19, 6, 6.0, 1.0, 0.0);//node 19        6.0    1.   0.
            theDomain.AddNode(node19);

            var node20 = new OpenSees.Components.NodeWrapper(20, 6, 6.0, 0.5, 0.0);//node 20        6.0   0.5   0.
            theDomain.AddNode(node20);

            var node21 = new OpenSees.Components.NodeWrapper(21, 6, 6.0, 0.0, 0.0);//node 21        6.0    0.   0.
            theDomain.AddNode(node21);

            var node22 = new OpenSees.Components.NodeWrapper(22, 6, 7.0, 1.0, 0.0);//node 22        7.0    1.   0.
            theDomain.AddNode(node22);

            var node23 = new OpenSees.Components.NodeWrapper(23, 6, 7.0, 0.5, 0.0);//node 23        7.0   0.5   0.
            theDomain.AddNode(node23);

            var node24 = new OpenSees.Components.NodeWrapper(24, 6, 7.0, 0.0, 0.0);//node 24        7.0    0.   0.
            theDomain.AddNode(node24);

            var node25 = new OpenSees.Components.NodeWrapper(25, 6, 8.0, 1.0, 0.0);//node 25        8.0    1.   0.
            theDomain.AddNode(node25);

            var node26 = new OpenSees.Components.NodeWrapper(26, 6, 8.0, 0.5, 0.0);//node 26        8.0   0.5   0.
            theDomain.AddNode(node26);

            var node27 = new OpenSees.Components.NodeWrapper(27, 6, 8.0, 0.0, 0.0);//node 27        8.0    0.   0.
            theDomain.AddNode(node27);

            var node28 = new OpenSees.Components.NodeWrapper(28, 6, 9.0, 1.0, 0.0);//node 28        9.0    1.   0.
            theDomain.AddNode(node28);

            var node29 = new OpenSees.Components.NodeWrapper(29, 6, 9.0, 0.5, 0.0);//node 29        9.0   0.5   0.
            theDomain.AddNode(node29);

            var node30 = new OpenSees.Components.NodeWrapper(30, 6, 9.0, 0.0, 0.0);//node 30        9.0    0.   0.
            theDomain.AddNode(node30);

            var node31 = new OpenSees.Components.NodeWrapper(31, 6, 10.0, 1.0, 0.0);//node 31         10.    1.   0.
            theDomain.AddNode(node31);

            var node32 = new OpenSees.Components.NodeWrapper(32, 6, 10.0, 0.5, 0.0);//node 32         10.   0.5   0.
            theDomain.AddNode(node32);

            var node33 = new OpenSees.Components.NodeWrapper(33, 6, 10.0, 0.0, 0.0);//node 33         10.    0.   0.
            theDomain.AddNode(node33);

            foreach (var tag in new[] { 1, 2, 3 })
                for (var dof = 0; dof < 6; dof++)
                    theDomain.AddSP_Constraint(new OpenSees.Components.Constraints.SP_ConstraintWrapper(tag, dof, 0.0, true));

            var elasticMat = new OpenSees.Materials.NDMaterials.ElasticIsotropicMaterialWrapper(2, 1.2e6, 0, 0);
            var plateFiberMat = new OpenSees.Materials.NDMaterials.PlateFiberMaterialWrapper(601, elasticMat);
            var section = new OpenSees.Materials.Sections.MembranePlateFiberSectionWrapper(701, 0.1, plateFiberMat);

            //element ShellNLDKGQ  1  2  5  4  1 701
            var ele1 = new OpenSees.Elements.ShellNLDKGQWrapper(1, 2, 5, 4, 1, section);
            theDomain.AddElement(ele1);

            //element ShellNLDKGQ  2  3  6  5  2 701
            var ele2 = new OpenSees.Elements.ShellNLDKGQWrapper(2, 3, 6, 5, 2, section);
            theDomain.AddElement(ele2);

            //element ShellNLDKGQ  3   5  8  7  4 701
            var ele3 = new OpenSees.Elements.ShellNLDKGQWrapper(3, 5, 8, 7, 4, section);
            theDomain.AddElement(ele3);

            //element ShellNLDKGQ  4   6  9  8  5 701
            var ele4 = new OpenSees.Elements.ShellNLDKGQWrapper(4, 6, 9, 8, 5, section);
            theDomain.AddElement(ele4);

            //element ShellNLDKGQ  5   8 11 10  7 701
            var ele5 = new OpenSees.Elements.ShellNLDKGQWrapper(5, 8, 11, 10, 7, section);
            theDomain.AddElement(ele5);

            //element ShellNLDKGQ  6   9 12 11  8 701
            var ele6 = new OpenSees.Elements.ShellNLDKGQWrapper(6, 9, 12, 11, 8, section);
            theDomain.AddElement(ele6);

            //element ShellNLDKGQ  7  11 14 13 10 701
            var ele7 = new OpenSees.Elements.ShellNLDKGQWrapper(7, 11, 14, 13, 10, section);
            theDomain.AddElement(ele7);

            //element ShellNLDKGQ  8  12 15 14 11 701
            var ele8 = new OpenSees.Elements.ShellNLDKGQWrapper(8, 12, 15, 14, 11, section);
            theDomain.AddElement(ele8);

            //element ShellNLDKGQ  9  14 17 16 13 701
            var ele9 = new OpenSees.Elements.ShellNLDKGQWrapper(9, 14, 17, 16, 13, section);
            theDomain.AddElement(ele9);

            //element ShellNLDKGQ  10  15 18 17 14 701
            var ele10 = new OpenSees.Elements.ShellNLDKGQWrapper(10, 15, 18, 17, 14, section);
            theDomain.AddElement(ele10);

            //element ShellNLDKGQ  11  17 20 19 16 701
            var ele11 = new OpenSees.Elements.ShellNLDKGQWrapper(11, 17, 20, 19, 16, section);
            theDomain.AddElement(ele11);

            //element ShellNLDKGQ  12  18 21 20 17 701
            var ele12 = new OpenSees.Elements.ShellNLDKGQWrapper(12, 18, 21, 20, 17, section);
            theDomain.AddElement(ele12);

            //element ShellNLDKGQ  13  20 23 22 19 701
            var ele13 = new OpenSees.Elements.ShellNLDKGQWrapper(13, 20, 23, 22, 19, section);
            theDomain.AddElement(ele13);

            //element ShellNLDKGQ  14  21 24 23 20 701
            var ele14 = new OpenSees.Elements.ShellNLDKGQWrapper(14, 21, 24, 23, 20, section);
            theDomain.AddElement(ele14);

            //element ShellNLDKGQ  15  23 26 25 22 701
            var ele15 = new OpenSees.Elements.ShellNLDKGQWrapper(15, 23, 26, 25, 22, section);
            theDomain.AddElement(ele15);

            //element ShellNLDKGQ  16  24 27 26 23 701
            var ele16 = new OpenSees.Elements.ShellNLDKGQWrapper(16, 24, 27, 26, 23, section);
            theDomain.AddElement(ele16);

            //element ShellNLDKGQ  17  26 29 28 25 701
            var ele17 = new OpenSees.Elements.ShellNLDKGQWrapper(17, 26, 29, 28, 25, section);
            theDomain.AddElement(ele17);

            //element ShellNLDKGQ  18  27 30 29 26 701
            var ele18 = new OpenSees.Elements.ShellNLDKGQWrapper(18, 27, 30, 29, 26, section);
            theDomain.AddElement(ele18);

            //element ShellNLDKGQ  19  29 32 31 28 701
            var ele19 = new OpenSees.Elements.ShellNLDKGQWrapper(19, 29, 32, 31, 28, section);
            theDomain.AddElement(ele19);

            //element ShellNLDKGQ  20  30 33 32 29 701
            var ele20 = new OpenSees.Elements.ShellNLDKGQWrapper(20, 30, 33, 32, 29, section);
            theDomain.AddElement(ele20);

            #region recorders
            var savepath = System.Environment.CurrentDirectory + @"\opsnet_results-exp7\";
            if (!System.IO.Directory.Exists(savepath))
                System.IO.Directory.CreateDirectory(savepath);

            {
                var nodeTags = new int[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33 };
                {
                    var opsstream = new OpenSees.Handlers.DataFileStreamWrapper(savepath + $@"\disp{1}.opsnet.txt");
                    var recorder = new OpenSees.Recorders.NodeRecorderWrapper(new OpenSees.IDWrapper(new int[] { 0, }), new OpenSees.IDWrapper(nodeTags), 0, "disp", theDomain, opsstream);
                    theDomain.AddRecorder(recorder);
                }

                {
                    var opsstream = new OpenSees.Handlers.DataFileStreamWrapper(savepath + $@"\disp{3}.opsnet.txt");
                    var recorder = new OpenSees.Recorders.NodeRecorderWrapper(new OpenSees.IDWrapper(new int[] { 2, }), new OpenSees.IDWrapper(nodeTags), 0, "disp", theDomain, opsstream);
                    theDomain.AddRecorder(recorder);
                }

                {
                    var opsstream = new OpenSees.Handlers.DataFileStreamWrapper(savepath + $@"\reaction{1}.opsnet.txt");
                    var recorder = new OpenSees.Recorders.NodeRecorderWrapper(new OpenSees.IDWrapper(new int[] { 0,1,2 }), new OpenSees.IDWrapper(nodeTags), 0, "reaction", theDomain, opsstream);
                    theDomain.AddRecorder(recorder);
                }
            }

            #endregion

            #region loading 
            var theSeries = new OpenSees.Components.Timeseries.LinearSeriesWrapper();
            var theLoadPattern = new OpenSees.Components.LoadPatterns.LoadPatternWrapper(1);
            theLoadPattern.SetTimeSeries(theSeries);
            theDomain.AddLoadPattern(theLoadPattern);

            var theLoadValues = new OpenSees.VectorWrapper(6);
            theLoadValues.Zero();
            theLoadValues[4] = 15.702963;

            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(1, 32, theLoadValues, false), 1);


            var theModel = new OpenSees.AnalysisModelWrapper();
            var theSolnAlgo = new OpenSees.Algorithms.KrylovNewtonWrapper();
            var theIntegrator = new OpenSees.Integrators.Static.LoadControlWrapper(0.001, 1, 0.001, 0.001);
            var theHandler = new OpenSees.Handlers.TransformationConstraintHandlerWrapper();
            var theRCM = new OpenSees.GraphNumberers.RCMWrapper(false);
            var theNumberer = new OpenSees.Numberers.DOF_NumbererWrapper(theRCM);
            var theSolver = new OpenSees.Systems.Linears.BandGenLinLapackSolverWrapper();
            var theSOE = new OpenSees.Systems.Linears.BandGenLinSOEWrapper(theSolver);
            var test = new OpenSees.ConvergenceTests.CTestNormDispIncrWrapper(1e-3, 1000, 2, 2, 1.0e10);
            var theAnalysis = new OpenSees.Analysis.StaticAnalysisWrapper(
                theDomain,
                theHandler,
                theNumberer,
                theModel,
                theSolnAlgo,
                theSOE,
                theIntegrator,
                test);

            theAnalysis.Analyze(4000);


            Console.ReadKey();
            #endregion
        }