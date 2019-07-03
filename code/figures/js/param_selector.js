// Define the constants;
var inducer = c_range;
var nns = 4600000;
var n = 2;
var operators = ['O1', 'O2', 'O3', 'Oid'];

// Define sliders and inputs, and 
var input_option = Radio.active; 
var drop_option = Drop.value;
var fit_option = Fit.value;
var description = Desc;
var fit_description = FitDesc;
var epri = 0;
var repressors = [22, 60, 124, 260, 1220, 1740];
var operators = ['O1', 'O2', 'O3', 'Oid']

// Case 1 -> Use numeric input
if  (input_option == 0 ) {
    var energies = [-15.3, -13.9, -9.7, -17]
    var ka_val = parseFloat(Ka_numer.value); 
    var ki_val = parseFloat(Ki_numer.value); 
    var krr_val = Math.exp(-parseFloat(Krr_numer.value)); 
    var epri = 0;
}

// Case I -> Choose literature values
else if (input_option == 1) {
    var energies = [-15.3, -13.9, -9.7, -17]
    if (drop_option == 'razo') {
        var ka_val = 139;
        var ki_val = 0.53;
        var krr_val = 0.01;
        description.text = "<b> Razo-Mejia <i>et al.</i>, Cell Systems, 2018</b><br/> Ka = 139 µM<br/> Ki = 0.53 µM<br/> Δε_AI= 4.5 kT<br/> Inactive repressor DNA binding is negligible";
        var epri = 0;
    }

    else if (drop_option == 'one_not_enough_x0') {
        var ka_val = 60;
        var ki_val = 4.2; 
        var krr_val = 2.3;
        description.text = "<b>Daber, Sharp, and Lewis. J. Mol. Biol. 2009</b><br/>Ka = 60 µM<br/> Ki = 4 µM<br/> ∆ε_AI = -0.83 kT<br/> Inactive repressor DNA binding is negligible"
        var epri = 0
    }

    else if (drop_option == 'one_not_enough_x100000') {
        var ka_val = 60;
        var ki_val = 4; 
        var krr_val = 2.3;
        var epri = Math.log(1 / 100000);
        description.text = "<b>Daber, Sharp, and Lewis. J. Mol. Biol. 2009</b><br/>Ka = 60 µM<br/> Ki = 4 µM<br/> Δε_AI = -0.83 kT <br/> Inactive repressor DNA binding 1 / 100,000 strength of active"
    }

    else if (drop_option == 'one_not_enough_x10000') {
        var ka_val = 60;
        var ki_val = 4; 
        var krr_val = 2.4;
        var epri = Math.log(1 / 10000);
        description.text = "<b>Daber, Sharp, and Lewis. J. Mol. Biol. 2009</b><br/>Ka = 60 µM<br/> Ki = 4 µM<br/> Δε_AI = -0.87 kT <br/> Inactive repressor DNA binding 1 / 10,000 strength of active"
    }

    else if (drop_option == 'one_not_enough_x1000') {
        var ka_val = 79;
        var ki_val = 4; 
        var krr_val = 1.3;
        var epri = Math.log(1 / 1000);
        description.text = "<b>Daber, Sharp, and Lewis. J. Mol. Biol. 2009</b><br/>Ka = 79 µM<br/> Ki = 4 µM<br/> ∆ε_AI = -0.87 kT <br/> Inactive repressor DNA binding 1 / 1000 strength of active"
    }

    else if (drop_option == 'one_not_enough_x100') {
        var ka_val = 812;
        var ki_val = 0.3; 
        var krr_val = 0.01;
        var epri = Math.log(1 / 100);
        description.text = "<b>Daber, Sharp, and Lewis. J. Mol. Biol. 2009</b><br/>Ka = 79 µM<br/> Ki = 4 µM<br/> Δε_AI = -4.6 kT <br/> Inactive repressor DNA binding 1 / 100 strength of active"
    }

    else if (drop_option == 'daber_muts')  {
        var ka_val =  16;
        var ki_val = 2;
        var krr_val = 5.8;
        var epri = 0
        description.text = "<b> Daber, Sochor, and  Lewis. J. Mol. Biol. 2011</b><br/>Ka = 16 µM<br/> Ki = 2 µM<br/> Δε_AI = -1.75 kT<br/> Inactive repressor DNA binding is negligible"
    }

    else if (drop_option == 'ogorman') {
        var ka_val = 133;
        var ki_val = 4;
        var krr_val = 0.7;
        var epri = Math.log(1 / 1000);
        description.text = "<b>O'Gorman <i> et al.</i>, JBC, 1980</b><br/>Ka = 133 µM<br/> Ki = 4 µM<br/> ∆ε_AI = 0.35 kT <br/> Inactive repressor DNA binding 1 / 1000 strength of active"
    }
}


// Case 2: Choose global fit values. 
else if (input_option==2) { 
    var epri = 0;

    if (fit_option == 'razo') { 
    var energies = [-15.121, -13.4, -9.21, -17.1]
    var ka_val = 225;
    var ki_val = 0.81
    var krr_val = Math.exp(-4.5);
    fit_description.text = "<b>Global fit using Razo-Mejia <i> et al.</i>, Cell Sys.,  2018 ∆ε_AI = 4.5 kT </b><br/>∆ε_RA (O1 operator) = -15.1 \u00B1 0.1 kT <br/> ∆ε_RA (O2 operator) = -13.4 \u00B1 0.1 kT<br/> ∆ε_RA (O3 operator) = -9.21 \u00B1 0.06 kT<br/> ∆ε_RA (Oid operator) = -17 \u00B1 5 kT<br/> Ka = 225 \u00B1 20 µM<br/> Ki = 0.81 \u00B1 0.05 µM<br/>"

    } 

    if (fit_option == 'ogorman') {
    var energies = [-15.7, -14.0, -9.85, -19];
    var ka_val = 267;
    var ki_val = 5.45;
    var krr_val = Math.exp(-0.35);
    fit_description.text = "<b>Global fit using O'Gorman <i> et al.</i>, JBC, 1980 ∆ε_AI = 0.35 kT </b><br/>∆ε_RA (O1 operator) = -15.7 \u00B1 0.1 kT <br/> ∆ε_RA (O2 operator) = -13.4 \u00B1 0.1 kT<br/> ∆ε_RA (O3 operator) = -9.85 \u00B1 0.06 kT<br/> ∆ε_RA (Oid operator) = -19 \u00B1 5 kT<br/> Ka = 270 \u00B1 20 µM<br/> Ki = 5.5 \u00B1 0.4 µM<br/>"
    }

    else if (fit_option == 'daber_muts') {
    var energies = [-17.111, -15.4, -11.29, -20]
    var ka_val = 288.92;
    var ki_val = 8.1927;
    var krr_val = Math.exp(1.75785); 
    fit_description.text = "<b>Global fit using Daber <i> et al.</i>, J. Mol. Biol., 2011 ∆ε_AI = -1.75 kT </b><br/>∆ε_RA (O1 operator) = -17.1 \u00B1 0.1 kT <br/> ∆ε_RA (O2 operator) = -13.4 \u00B1 0.1 kT<br/> ∆ε_RA (O3 operator) = -9.85 \u00B1 0.06 kT<br/> ∆ε_RA (Oid operator) = -20 \u00B1 5 kT<br/> Ka = 270 \u00B1 20 µM<br/> Ki = 8.2 \u00B1 0.4 µM<br/>"
    }
}

// Define the function for computing the fold-change
function foldChange(R, epr, epri) {
    fc = []
    if (epri != 0) {
        var pref = 1;
    }
    else {
        var pref = 0
    }
    for (var i = 0; i < inducer.length; i++) {
    var numer = Math.pow(1 + inducer[i] / ka_val, n);
    var denom = numer + krr_val * Math.pow(1 + inducer[i] / ki_val, n);
    fc[i] =  1 / (1 + (numer / denom) * (R / nns) * Math.exp(-epr) + pref * (1 - (numer / denom)) * (R / nns) * Math.exp(-(epr - epri)));
    }
    return fc;
}

function updateFoldChange(source, epr) {
    var ys = [];
    var cs = [];
    for (var i = 0; i < repressors.length; i++) {
        ys[i] = foldChange(repressors[i], epr, epri);
        cs[i] = repressors[i];
    }
    source.data['fold_change'] = ys;
    source.change.emit()
}

function updateLeakiness(source) {
    var ys=[]; 
    var pact = 1 / (1 + krr_val);
        if (epri != 0) {
            var pref = 1;
        }
        else {
            var pref = 0;
        }
    for (var i =0 ; i < operators.length; i++) {
        var _ys = []
        for (var j=0; j < rep_range.length; j++) {
            var ra = rep_range[j] * pact / nns;
            var ri = rep_range[j] * (1 - pact) / nns;

            _ys[j] = 1 / (1 + ra * Math.exp(-energies[i]) + pref * ri * Math.exp(-(energies[i] - epri)));
            console.log(_ys[j])
    }
    ys[i] = _ys;
    }
    // source.data['repressors'] = xs;
    source.data['leakiness'] = ys;
    source.change.emit();
}


// Define some functions for computing quantities.
function computeFugacity(ep,  epri, M) {
    // Compute the residual active probability. 
    var pact = (1 / (1 + krr_val));
    var fug = [];
    if (epri != 0) {
        var pref = 1;
    }

    else {
        var pref = 0;
    }
    var x = Math.exp(-ep);
    var xi = Math.exp(-(ep - epri));
    for (var j = 0; j < rep_range.length; j++) {
        var ra = rep_range[j] * pact;
        var ri = rep_range[j] * (1 - pact);
        var B_act = (ra * (1 + x) - M * x - 4600000) 
        var B_inact = (ri * (1 + xi) - M * xi - 4600000) 
        var C_act = ra;
        var C_inact = ri;
        var A_act = x * (ra - M - 4600000);
        var A_inact = xi * (ri - M - 4600000);
        var numer_act = -B_act - Math.sqrt(Math.pow(B_act, 2) - 4 * A_act * C_act);
        var numer_inact = -B_inact - Math.sqrt(Math.pow(B_inact, 2) - 4 * A_inact * C_inact);
        var denom_act = 2 * A_act;
        var denom_inact = 2 * A_inact;
        fug[j] = 1 / (1 + (numer_act / denom_act) * x + pref * (numer_inact /denom_inact) * xi) 
 }
 return fug;
}

var M = [64, 52]
function updateFoldChangeFugacity(source) {
    ys = [];
    for (var j = 0; j < M.length; j++) {
        ys[j] = computeFugacity(energies[0], epri, M[j]);
    }
    source.data['leakiness'] = ys;
    source.change.emit();
    }
// Define the three sources
updateFoldChange(o1_source, energies[0], epri);
updateFoldChange(o2_source, energies[1], epri);
updateFoldChange(o3_source, energies[2], epri);
updateLeakiness(leak_source);
updateFoldChangeFugacity(fug_source);
