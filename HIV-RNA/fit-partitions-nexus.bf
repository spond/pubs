RequireVersion("2.26");


VERBOSITY_LEVEL = 0;
//LF_SMOOTHING_SCALER = 0.1;

// namespace 'utility' for convenience functions
LoadFunctionLibrary("libv3/UtilityFunctions.bf");


// namespace 'alignments' for alignment-related functions
LoadFunctionLibrary("libv3/tasks/alignments.bf");

// namespace 'trees' for trees-related functions
LoadFunctionLibrary("libv3/tasks/trees.bf");


// namespace 'io' for interactive/datamonkey i/o functions
LoadFunctionLibrary("libv3/IOFunctions.bf");

LoadFunctionLibrary("libv3/models/DNA/HKY85.bf");

// namespace 'estimators' for various estimator related functions
LoadFunctionLibrary("libv3/tasks/estimators.bf");

// namespace 'ancestral' for ancestral reconstruction functions
LoadFunctionLibrary("libv3/tasks/ancestral.bf");

// namespace 'stats' for various descriptive stats functions
LoadFunctionLibrary("libv3/stats.bf");


LoadFunctionLibrary("libv3/terms-json.bf");


/*------------------------------------------------------------------------------*/

utility.ToggleEnvVariable ("NORMALIZE_SEQUENCE_NAMES", 1);

io.DisplayAnalysisBanner({
    "info": "Use a partition annotation file, and a time-scaled MCC tree, split
    the alignment into partitions, fit partition-specific models with substitution rates
    partitioned into terminal and interior nodes, report the MLEs",
    "version": "1.00",
    "reference": "'_Combining chemical and evolutionary information to effectively evaluate HIV and SIV genetic heterogeneity and RNA structural conservation within and among infected hosts_', by Rife _et. al._",
    "authors": "Sergei L Kosakovsky Pond",
    "contact": "spond@temple.edu"
});


SetDialogPrompt ("Please load the NEXUS file with CHARSETs");

shape.file_info = alignments.ReadNucleotideDataSet ("shape.nuc_data", None);

DataSetFilter shape.all = CreateFilter (shape.nuc_data, 1);

if (None != shape.file_info["name-mapping"]) {
    shape.tree = trees.LoadAnnotatedTopologyAndMap(1,shape.file_info["name-mapping"]);
} else {
    shape.tree = trees.LoadAnnotatedTopology(1);
}


io.CheckAssertion ("None!=shape.file_info[\"partitions\"]", "The input file did not contain any partition information");
shape.partition = shape.file_info["partitions"];

io.ReportProgressMessageMD ("", "data", "Loaded an MSA with **" + shape.file_info ["sequences"] + "** sequences and **" + shape.file_info ["sites"] + "** sites from \`" + shape.file_info["file"] + "\`");
io.ReportProgressMessageMD ("", "data" "Read **" + Abs (shape.partition) + "** partitions on **" + shape.file_info["sites"] + "** sites.");

partition_fits = {};

shape.model_assignment = {"leaves": {}, "inodes" : {}};

function do_partition_tree (key, value) {
    if (value == "leaf") {
        (shape.model_assignment["leaves"]) [key] = 1;
    } else {
        (shape.model_assignment["inodes"]) [key] = 1;
    }
}

(shape.tree["partitioned"])["do_partition_tree"][""];
utility.ToggleEnvVariable("USE_LAST_RESULTS", 1);

logL_by_partition = {};

for (partition_id = 0; partition_id < Rows (shape.partition); partition_id += 1) {

    parameters.DeclareGlobal ("leaf_rate", None);
    parameters.DeclareGlobal ("inode_rate", None);
    partition_name = shape.partition [partition_id][0];
    DataSetFilter partition_filter = CreateFilter (shape.nuc_data, 1, shape.partition [partition_id][1]);
    partition_length = partition_filter.sites;
    io.ReportProgressMessageMD ("", ""+partition_name, "Fitting partition " + (partition_name) + " with " + describe_partition ("partition_filter"));

    model_leaves = model.generic.DefineModel ("models.DNA.HKY85.ModelDescription", "model.leaves", {
        "0": "terms.global"
    }, "partition_filter", None);

    ((model_leaves["parameters"])[terms.global])["Leaf rate scaler"] = "leaf_rate";

    model_ibranches = model.generic.DefineModel ("models.DNA.HKY85.ModelDescription", "model.inodes", {
        "0": "terms.global"
    }, "partition_filter", None);


    ((model_ibranches["parameters"])[terms.global])["Internal branch rate scaler"] = "inode_rate";

    model.ApplyModelToTree("partition_tree", shape.tree, {
            "leaves": model_leaves ["id"], "inodes": model_ibranches ["id"]} , shape.model_assignment);


    // rename parameters in the leaf model

    parameters.SetRange ("leaf_rate", {"`terms.lower_bound`" : 0, "`terms.upper_bound`" : 1e-3});
    parameters.SetRange ("inode_rate", {"`terms.lower_bound`" : 0, "`terms.upper_bound`" : 1e-3});
    LikelihoodFunction partition_lf = (partition_filter, partition_tree);

    node_names = utility.Keys (shape.tree["partitioned"]);

    chrono_time = 0;

    utility.ForEachPair (shape.tree["partitioned"], "_node_", "_node_class_",
        'if (_node_class_ == "leaf") {
            chrono_time += (shape.tree[terms.json.attribute.branch_length])[_node_];
            estimators.applyBranchLength("partition_tree", _node_, model_leaves, {"`terms.branch_length_scaler`" : "leaf_rate", "`terms.branch_length`" : (shape.tree[terms.json.attribute.branch_length])[_node_]});
        } else {
            estimators.applyBranchLength("partition_tree", _node_, model_ibranches, {"`terms.branch_length_scaler`" : "inode_rate", "`terms.branch_length`" : (shape.tree[terms.json.attribute.branch_length])[_node_]});
        }
    ');



    Optimize (mles, partition_lf);
    results = estimators.ExtractMLEs("partition_lf", {
        "model.leaves": model_leaves,
        "model.inodes": model_ibranches
    });

    evol_time = 0;
    ((shape.model_assignment["inodes"]))["extract_inode_lengths"][""];



    //fprintf (stdout, inode_rate, ":", evol_time, ":", chrono_time, "\n", mles, "\n");
    io.ReportProgressMessageMD ("", ""+partition_name, "* log (L)  **" + format_number(mles[1][0]) + "**");
    io.ReportProgressMessageMD ("", ""+partition_name, "* &pi; (A)  **" + report_mle ((model_leaves[terms.efv_estimate])[0], None, None) + "**");
    io.ReportProgressMessageMD ("", ""+partition_name, "* &pi; (C)  **" + report_mle ((model_leaves[terms.efv_estimate])[1], None, None) + "**");
    io.ReportProgressMessageMD ("", ""+partition_name, "* &pi; (G)  **" + report_mle ((model_leaves[terms.efv_estimate])[2], None, None) + "**");
    io.ReportProgressMessageMD ("", ""+partition_name, "* &pi; (T)  **" + report_mle ((model_leaves[terms.efv_estimate])[3], None, None) + "**");
    io.ReportProgressMessageMD ("", ""+partition_name, "* Internal node rate  **" + report_mle (evol_time / chrono_time , parameters.GetProfileCI ("inode_rate", "partition_lf", 0.95), 1) + "**");
    io.ReportProgressMessageMD ("", ""+partition_name, "* Internal node kappa **" + report_mle (((results["global"])[terms.transition_transversion_ratio])["MLE"], parameters.GetProfileCI (((results["global"])[terms.transition_transversion_ratio])["ID"], "partition_lf", 0.95), 0) + "**");

    DeleteObject (partition_lf);
}

function extract_inode_lengths (nn, value) {
    bl = (((results[terms.json.attribute.branch_length])[0]) [nn])["MLE"];
    if (None != bl) {
        evol_time += bl;
    }
}

lfunction format_number (value) {
    return Format (value, 0, Max (3, -Log (Abs(value))/Log (10) $ 1 + 3));
}

function report_mle (value, ci, multiply) {
    if (None != ci) {
        if (multiply) {
            return "" + format_number(value) + " [" + format_number(value*ci[terms.lower_bound]/ci[terms.MLE]) + " - " + format_number(value*ci[terms.upper_bound]/ci[terms.MLE]) + "]";
        }
        return "" + format_number(value) + " [" +  format_number(ci[terms.lower_bound]) + " - " + format_number(ci[terms.upper_bound]) + "]";
    }
    return format_number(value);
}

lfunction describe_partition (filter_name) {

    GetInformation (info, ^filter_name);
    site_count           = Abs (info[0]);
    sequence_count       = Columns (info);
    constant             = 0;
    variable             = 0;
    phylogenetically_inf = 0;

    for (k = 0; k < site_count; k += 1) {
        site_profile = {};
        for (s = 0; s < sequence_count; s+=1) {
            site_profile[(info[s])[k]] += 1;
        }
        if (Abs (site_profile) == 1) {
            constant += 1;
        } else {
            variable += 1;
            pi = 0;
            chars = Rows (site_profile);
            for (c = 0; c < Abs (site_profile); c+=1) {
                pi += (site_profile[chars[c]]) > 1;
            }
            phylogenetically_inf += (pi > 1);
        }
    }
    return "**" + constant + "** constant sites, **" + variable + "** variable sites, and **" +  phylogenetically_inf + "** phylogenetically informative sites";
}
