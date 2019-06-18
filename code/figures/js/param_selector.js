// Define data sources
var o1 = o1_source.data;
var o2 = o2_source.data;
var o3 = o3_source.data;

// Define the constants;
var inducer = c_range;
var nns = 4600000;
var n = 2;
var repressors = [22, 60, 124, 260, 1220, 1740];

// Define sliders and inputs, and 
var input_option = Radio.active; 

if  (input_option == 0 ) {
    var ka_val = parseFloat(Ka_numer.value); 
    var ki_val = parseFloat(Ki_numer.value); 
    var krr_val = parseFloat(Krr_numer.value); 

}

else {
    var ka_val = Ka_slider.value; 
    var ki_val = Ki_slider.value;
    var krr_val = Krr_slider.value; 
}

// Define the function for computing the fold-change
function foldChange(R, epr) {
    fc = []
    for (var i = 0; i < inducer.length; i++) {
    var numer = Math.pow(1 + inducer[i] / ka_val, n);
    var denom = numer + krr_val * Math.pow(1 + inducer[i] / ki_val, n);
    fc[i] =  1 / ( 1 + (numer / denom) * (R / nns) * Math.exp(-epr));
    }
    return fc;
}

function updateFoldChange(source, epr) {
    var xs = [];
    var ys = [];
    var cs = []
    for (var i = 0; i < repressors.length; i++) {
        xs[i] = inducer;
        ys[i] = foldChange(repressors[i], epr);
        cs[i] = repressors[i];
    }
    console.log(ys)
    source.data['IPTGuM'] = xs;
    source.data['fold_change'] = ys;
    return source
}

// Define the three sources
o1_source = updateFoldChange(o1_source, -15.3);
o2_source = updateFoldChange(o2_source, -13.9);
o3_source = updateFoldChange(o3_source, -9.7);
o1_source.change.emit();
o2_source.change.emit();
o3_source.change.emit();
