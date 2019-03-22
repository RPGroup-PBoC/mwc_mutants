function linspace(start, stop, n) {
    var x = [];
    var currValue = start;
    var step = (stop - start) / (n - 1);
    for (var i = 0; i < n; i++) {
        x.push(currValue);
        currValue += step;
    }
    return x;
}

function logspace(start, stop, n) {
    var x = [];
    var currValue = start;
    var step = (stop - start) / n;
    for (var i = 0; i < n; i++) {
        x.push(Math.pow(10, currValue));
        currValue += step   
    }
    return x
}

function derivEpAI(epValue, c, ka, ki, n){
    var allo_ratio =  Math.pow((ka * (ki + c)) / (ki * (ka + c)), n);
    return allo_ratio / (allo_ratio + Math.exp(epValue));
}

function derivC(ep, cValue,  ka, ki, n){
    var numer = -n * Math.pow((ka * (ki + cValue))/(ki * (ka + cValue)),n) * (ka - ki);
    var denom_a = (ka + cValue) * (ki + cValue);
    var denom_b = (Math.pow((ka * (ki + cValue))/ (ki * (ka + cValue)), n) + Math.exp(ep));
    return numer / (denom_a * denom_b);
}

function derivKa(ep, c, kaValue, ki, n) {
    var numer = -c * n * Math.pow((kaValue * (ki + c))/ (ki * (kaValue + c)), n);
    var denom_a = kaValue * (kaValue + c)
    var denom_b = Math.pow((kaValue * (ki + c))/ (ki * (kaValue + c)), n) + Math.exp(ep);
    return numer / (denom_a * denom_b);
}

function derivKi(ep, c, ka, kiValue, n) { 
    var numer = c * n * Math.pow((ka * (kiValue + c))/ (kiValue * (ka + c)), n);
    var denom_a = kiValue * (kiValue + c)
    var denom_b = Math.pow((ka * (kiValue + c)) / (kiValue * (ka + c)), n)  + Math.exp(ep)
    return numer / (denom_a * denom_b);
}

// Define the parameter values
var data = source.data;
var ka = Math.pow(10, Ka.value);
var ki = Math.pow(10, Ki.value);
var epAI = EpAI.value;
var c = Math.pow(10, C.value);
var n = nSites;

// Define the observed ranges
data['epAI'] = linspace(xepAI.start, xepAI.end, n_points); 
data['ka'] = logspace(Math.log10(xKa.start), Math.log10(xKa.end), n_points);
data['ki'] = logspace(Math.log10(xKi.start), Math.log10(xKi.end), n_points);
data['c']= logspace(Math.log10(xC.start), Math.log10(xC.end), n_points);
console.log(data['ka'])
// Call the functions to compute the derivatives
for (var i = 0; i < n_points; i++) {
   data['dF_depAI'][i] = derivEpAI(data['epAI'][i], c, ka, ki, n);
   data['dF_dc'][i] = derivC(epAI, data['c'][i], ka, ki, n);
   data['dF_dka'][i] = derivKa(epAI, c, data['ka'][i], ki, n);
   data['dF_dki'][i] = derivKi(epAI, c, ka, data['ki'][i], n);}

// Update the sources. 
source.change.emit()
