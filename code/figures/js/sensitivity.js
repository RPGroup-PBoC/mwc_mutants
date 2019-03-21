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
}

function derivEpAI(epValue, c, ka, ki, n){
    var allo_ratio =  Math.pow(1 + c / ka, n) / Math.pow(1 + c/ki, n);
    return Math.pow(1 + Math.exp(epValue) * allo_ratio, -1);
}

function derivC(ep, cValue,  ka, ki, n){
    var numer = n * Math.pow((ki * (ka + cValue))/(ka * (ki + cValue)),n) * (ka - ki);
    var denom_a = (ka + cValue) * (ki + cValue);
    var denom_b = (Math.pow((ki * (ka + cValue) / (ka * (ki + cValue))), n) + Math.exp(ep));
    return numer / (denom_a * denom_b);
}

function derivKa(ep, c, kaValue, ki, n) {
    var numer = c * n * Math.pow((ki * (kaValue + c))/ (kaValue * (ki + c)), n);
    var denom_a = kaValue * (kaValue + c)
    var denom_b = Math.pow((ki * (ka + c))/ (ka * (ki + c)), n) + Math.exp(ep);
    return numer / (denom_a * denom_b);
}

function derivKi(ep, c, ka, kiValue, n) { 
    var numer = -c * n * Math.pow((kiValue * (ka + c))/ (ka * (kiValue + c)), n);
    var denom_a = kiValue * (kiValue + c)
    var denom_b = Math.pow((kiValue * (ka + c)) / (ka * (kiValue + c)), n)  + Math.exp(ep)
    return numer / (denom_a * denom_b);
}

// Define the parameter values
var data = source.data;
var ka = ka.value;
var ki = ki.value;
var epAI = epAI.value;
var c = c.value;
var n = nSites.value;

// Enforce positivity for log-scaled variables. 
if (xKa.start <= 0) {
    xKa.start = 0.00000002
}
if (xKa.end <= 0) {
    xKa.end = 0.00000002
}
if (xKi.start <= 0) {
    xKi.start = 0.00000002
}
if (xKi.end <= 0) {
    xKi.end = 0.00000002
}
if (xC.start <= 0) {
    xC.start = 0.00000002
}
if (xC.end <= 0) {
    xC.end = 0.00000002
}
// Define the observed ranges
data['epAI'] = linspace(xepAI.start, xepAI.end, n_points); 
data['ka'] = linspace(xKa.start, xKa.end, n_points);
data['ki'] = linspace(xKi.start, xKi.end, n_points);
data['c'] = linspace(xC.start, xC.end, n_points);


// Call the functions to compute the derivatives
for (var i = 0; i < n_points; i++) {
   data['dF_depAI'][i] = derivEpAI(data['epAI'][i], c, ka, ki, n);
   data['dF_dc'][i] = derivC(epAI, data['c'][i], ka, ki, n);
   data['dF_dka'][i] = derivKa(epAI, c, data['ka'][i], ki, n);
   data['dF_dki'][i] = derivKi(epAI, c, ka, data['ki'][i], n);}

// Update the sources. 
source.change.emit()
