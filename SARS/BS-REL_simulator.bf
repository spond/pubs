CAP_OMEGA3_AT_THIS_VALUE = 10;
CAP_T1_AT_THIS_VALUE     = 10000;

LoadFunctionLibrary ("GrabBag");

SetDialogPrompt ("Please load a branch-site REL model to simulate data from");
ExecuteAFile    (PROMPT_FOR_FILE);
//ExecuteAFile    ("rdrp.bsrel.omega3.results.fit");
//ExecuteAFile ("100.0.fit.fit");

//capValues ("mixtureTree", CAP_OMEGA3_AT_THIS_VALUE, CAP_T1_AT_THIS_VALUE);

fprintf (stdout, "Extracting fitted tree information\n");

_mleLengths   = reportTreeLength   ("mixtureTree", 1);
desired_scale = prompt_for_a_value ("Desired total tree length:" ,10,0.001,1e10,0);

scaleTree ("mixtureTree", +((_mleLengths["Props"])[-1][1]), desired_scale, {"0": "t1"});
reportTreeLength ("mixtureTree", 1);

replicates = prompt_for_a_value ("How many replicates should be generated:" ,100,1,100000,1);

SetDialogPrompt ("Save replicates to this file.xx");
fprintf         (PROMPT_FOR_FILE, CLEAR_FILE, KEEP_OPEN, "Replicate,Path");
_outPrefix = LAST_FILE_PATH;

for (_i = 0; _i < replicates; _i += 1) {
    DataSet       sim           = SimulateDataSet (three_LF);
    DataSetFilter simFilter     = CreateFilter (sim,1);
    _saveTo = _outPrefix + "." + _i;
    fprintf (_outPrefix, "\n", _i, ",", _saveTo);
    fprintf (_saveTo, simFilter);
    SetParameter (STATUS_BAR_STATUS_STRING, "Replicate " + (_i+1) + " of " + replicates + " generated", 0);
}

fprintf (_outPrefix, CLOSE_FILE);

//--------------------------------------------------------------------------------------------------------------

function reportTreeLength (treeID, doPrint) {
    _mleLengths  = extractBranchInformation (treeID, "omega[0-9]+", "MGMatrix", "Paux", "codon3x4", 0);
    _totalLength = +((_mleLengths["Props"])[-1][1]);
    if (doPrint) {
        for (_i = 0; _i < Rows (_mleLengths["Props"]); _i += 1) {
            fprintf (stdout, (_mleLengths["Names"])[_i], "(", (_mleLengths["Props"])[_i][1],") represents ", Format((_mleLengths["Props"])[_i][1]/_totalLength * 100.,6,4), "% of the tree\n");
        }
        fprintf (stdout, "\nTotal fitted tree length = ", _totalLength , "\n");
    }
    return _mleLengths;
}

//--------------------------------------------------------------------------------------------------------------

function capValues (treeID, capO3, capT1) {
    _branchNames   = Eval ("BranchName(`treeID`,-1)");
    _branchCount   = Columns (_branchNames) - 1;
    for (_i = 0; _i < _branchCount; _i += 1) {
        ExecuteCommands (treeID + "." + _branchNames[_i] + ".omega3:<" + capO3);
        ExecuteCommands (treeID + "." + _branchNames[_i] + ".t1:<" + capT1);
    }
    return 0;
}


//--------------------------------------------------------------------------------------------------------------

function scaleTree (treeID, existingLength, desiredLength, parametersToScale) {
    _branchNames            = Eval ("BranchName(`treeID`,-1)");
    _multiplicativeFactor   = desiredLength/existingLength;
    _branchCount            = Columns (_branchNames) - 1;
    
    for (_i = 0; _i < _branchCount; _i += 1) {
        for (_p = 0; _p < Abs (parametersToScale); _p += 1) {
            ExecuteCommands (treeID + "." + _branchNames[_i] + "." + parametersToScale[_p] + "=_multiplicativeFactor*" + treeID + "." + _branchNames[_i] + "." + parametersToScale[_p]);
        }  
    }
    return None;
}

//--------------------------------------------------------------------------------------------------------------

function extractBranchInformation (treeID, _omega_pattern, _model_prefix, _freq_prefix, _freq_vector, _do_mult) {
    _branchNames   = Eval ("BranchName(`treeID`,-1)");
    _branchCount   = Columns (_branchNames) - 1;

    _toReturn         = {"Names": {}, "Props": {_branchCount,2}};  // cat count, branch length
    _brLenExpressions = {};
    _brLenLocals      = {};
    
    for (_i = 0; _i < _branchCount; _i += 1) {
        _expression = treeID + "\\." + _branchNames[_i] + "\\." + _omega_pattern;
        ExecuteCommands("GetInformation (_matchedVars, \"`_expression`\")");
        _catCount = Columns (_matchedVars);
        _toReturn ["Names"] + _branchNames[_i];
        (_toReturn["Props"])[_i][0] = _catCount;
        
        _brLen     = 0;
        _totalWeight = 0;
        fprintf (stdout, "\n", _branchNames[_i], "\n");
        for (_c = 1; _c <= _catCount; _c += 1) {
            _modelID = _model_prefix + _c;
            if (Abs (_brLenExpressions[_modelID]) == 0) {
                ExecuteCommands ("Model _internal = (`_modelID`, `_freq_vector`, _do_mult);");
                GetString (_bl, _internal, -1);
                _brLenExpressions [_modelID] = _bl;
                _localParameters = {};
                _p = 1;
                GetString (_locP, _internal, 0);
                while (_locP != None) {
                    _localParameters + _locP; 
                    GetString (_locP, _internal, _p);
                    _p += 1;
                }
                _brLenLocals [_modelID] = _localParameters;
            }
            for (_lpc = 0; _lpc < Abs (_brLenLocals [_modelID]); _lpc += 1) {
                ExecuteCommands ((_brLenLocals [_modelID])[_lpc] + "=" + treeID + "." + _branchNames[_i] + "." + (_brLenLocals [_modelID])[_lpc]);
            }
            _thisModelLength = Eval (_brLenExpressions[_modelID]);
            if (_c < _catCount) {
                _freqW = Eval (treeID + "." + _branchNames[_i] + "." + _freq_prefix + _c);
            } else {
                _freqW = 1-Eval (treeID + "." + _branchNames[_i] + "." + _freq_prefix + (_c-1));
            }
            _thisModelWeight = (1-_totalWeight)*_freqW;
            fprintf (stdout, "\n\tomega = ", Eval ("omega" + _c), ", weight ", _thisModelWeight);
            
            _totalWeight += _thisModelWeight;
            _brLen += _thisModelWeight * _thisModelLength;            
        }
        (_toReturn["Props"]) [_i][1] = _brLen/3;
    }
    
    
    return _toReturn;   
}