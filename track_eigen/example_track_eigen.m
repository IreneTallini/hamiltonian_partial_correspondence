p = ["cat0", "cat1", "cat2", "cat3", "cat4", "cat5", "cat6", "cat7", "cat8", "cat9"];
Ms = load_sequence_mesh(p, "TOSCA");
Ms = eigdec_multiple_meshes(Ms, 50);
track_eigen(Ms, 1:50);