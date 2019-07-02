// Define the constants;
var inducer = c_range;
var nns = 4600000;
var n = 2;
var operators = ['O1', 'O2', 'O3', 'Oid']
var energies = [-15.3, -13.9, -9.7, -17]
var repressors = [22, 60, 124, 260, 1220, 1740];

// Define sliders and inputs, and 
var input_option = Radio.active; 
var drop_option = Drop.value;
var description = Desc;
var epri = 0;

if  (input_option == 0 ) {
    var ka_val = parseFloat(Ka_numer.value); 
    var ki_val = parseFloat(Ki_numer.value); 
    var krr_val = parseFloat(Krr_numer.value); 
    var epri = 0

}

else if (input_option == 1) {
    var ka_val = Ka_slider.value; 
    var ki_val = Ki_slider.value;
    var krr_val = Krr_slider.value; 
    var epri = 0
}

else {
    if (drop_option == 'razo') {
        var ka_val = 139;
        var ki_val = 0.53;
        var krr_val = 0.01
        description.text = "<b> Razo-Mejia <i>et al.</i>, Cell Systems, 2018</b><br/> Ka = 139 µM<br/> Ki = 0.53 µM<br/> K_RR*=0.01<br/> Inactive repressor DNA binding is negligible";
        var epri = 0;
    }

    else if (drop_option == 'one_not_enough_x0') {
        var ka_val = 60;
        var ki_val = 4.2; 
        var krr_val = 2.3;
        description.text = "<b>Daber, Sharp, and Lewis. J. Mol. Biol. 2009</b><br/>Ka = 60 µM<br/> Ki = 4 µM<br/> K_RR* = 2.3<br/> Inactive repressor DNA binding is negligible"
        var epri = 0
    }

    else if (drop_option == 'one_not_enough_x100000') {
        var ka_val = 60;
        var ki_val = 4; 
        var krr_val = 2.3;
        var epri = Math.log(1 / 100000);
        description.text = "<b>Daber, Sharp, and Lewis. J. Mol. Biol. 2009</b><br/>Ka = 60 µM<br/> Ki = 4 µM<br/> K_RR* = 2.3<br/> Inactive repressor DNA binding 1 / 100,000 strength of active"
    }

    else if (drop_option == 'one_not_enough_x10000') {
        var ka_val = 60;
        var ki_val = 4; 
        var krr_val = 2.4;
        var epri = Math.log(1 / 10000);
        description.text = "<b>Daber, Sharp, and Lewis. J. Mol. Biol. 2009</b><br/>Ka = 60 µM<br/> Ki = 4 µM<br/> K_RR* = 2.4<br/> Inactive repressor DNA binding 1 / 10,000 strength of active"
    }

    else if (drop_option == 'one_not_enough_x1000') {
        var ka_val = 79;
        var ki_val = 4; 
        var krr_val = 1.3;
        var epri = Math.log(1 / 1000);
        description.text = "<b>Daber, Sharp, and Lewis. J. Mol. Biol. 2009</b><br/>Ka = 79 µM<br/> Ki = 4 µM<br/> K_RR* = 2.4<br/> Inactive repressor DNA binding 1 / 1000 strength of active"
    }

    else if (drop_option == 'one_not_enough_x100') {
        var ka_val = 812;
        var ki_val = 0.3; 
        var krr_val = 0.01;
        var epri = Math.log(1 / 100);
        description.text = "<b>Daber, Sharp, and Lewis. J. Mol. Biol. 2009</b><br/>Ka = 79 µM<br/> Ki = 4 µM<br/> K_RR* = 0.01<br/> Inactive repressor DNA binding 1 / 100 strength of active"
    }

    else if (drop_option == 'daber_muts')  {
        var ka_val =  16;
        var ki_val = 2;
        var krr_val = 5.8;
        var epri = 0
        description.text = "<b> Daber, Sochor, and  Lewis. J. Mol. Biol. 2011</b><br/>Ka = 16 µM<br/> Ki = 2 µM<br/> K_RR* = 5.8<br/> Inactive repressor DNA binding is negligible"
    }

    else if (drop_option == 'ogorman') {
        var ka_val = 133;
        var ki_val = 4;
        var krr_val = 0.7;
        var epri = Math.log(1 / 1000);
        description.text = "<b>O'Gorman <i> et al.</i>, JBC, 1980</b><br/>Ka = 133 µM<br/> Ki = 4 µM<br/> K_RR* = 0.7<br/> Inactive repressor DNA binding 1 / 1000 strength of active"
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

            _ys[j] = 1 / (1 + ra * Math.exp(-energies[i]) + pref * ri * Math.exp(-epri * energies[i]));
            console.log(_ys[j])
    }
    ys[i] = _ys;
    }
    // source.data['repressors'] = xs;
    source.data['leakiness'] = ys;
    source.change.emit();
}


// Define function to compute the fugacity. 
function computePact() {
    var numer = Math.pow(1 + 0 / ki_val, n);
    var denom = Math.pow(1 + 0 / ka_val, n);
    var prob = Math.pow(1 + krr_val* numer / denom, -1);
    return prob;
    }

function computeFugacity(pact, ep,  epri, M) {
    var fug = [];
    if (epri != 0) {
        var pref = 1;
    }

    else {
        var pref = 0;
    }
    var x = Math.exp(15.3);
    var xi = Math.exp(-epri * (-15.3));
    for (var j = 0; j < rep_range.length; j++) {

        var ra = rep_range[j] * pact;
        var ri = rep_range[j] * (1 - pact);
        var B_act = (ra * (1 + x) - M * x - nns) 
        var B_inact = (ri * (1 + xi) - M * xi - nns) 
        var C_act = ra;
        var C_inact = ri;
        var A_act = x * (ra - M - nns);
        var A_inact = xi * (ri - M - nns);
        var numer_act = -B_act - Math.sqrt(B_act*B_act - 4 * A_act * C_act);
        var numer_inact = -B_inact - Math.sqrt(B_inact*B_inact - 4 * A_inact * C_inact);
        var denom_act = 2 * A_act;
        var denom_inact = 2 * A_inact;
        fug[j] = 1 / (1 + (numer_act / denom_act) * x + pref * (numer_inact /denom_inact) * xi) 
 }
 return fug;
}

var M = [64, 52]
function updateFoldChangeFugacity(source) {
    var pact = 1 / (1 + krr_val);
    ys = [];
    for (var j = 0; j < M.length; j++) {
        ys[j] = computeFugacity(pact, -15.3, epri, M[j]);
    }
    source.data['leakiness'] = ys;
    source.change.emit();
    }
// Define the three sources
updateFoldChange(o1_source, -15.3, epri);
updateFoldChange(o2_source, -13.9, epri);
updateFoldChange(o3_source, -9.7, epri);
updateLeakiness(leak_source);
updateFoldChangeFugacity(fug_source);
