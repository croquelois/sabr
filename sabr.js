(function(undefined){
var globalScope = (typeof global !== 'undefined') ? global : this;
var hasModule = (typeof module !== 'undefined') && module.exports;
if(hasModule) numeric = require('./numeric');

function sabr(p,isNormal){
  var ln = Math.log;
  var sqrt = Math.sqrt;
  var pow = Math.pow;
  var abs = Math.abs;
  var exp = Math.exp;
  var rndUni = Math.random;
  
  function sq(v){ return v*v; }
  function D(theta){
    return ln((sqrt(1-2*p.rho*theta+theta*theta)+theta-p.rho)/(1-p.rho));
  }
  function gamma(K){
    var beta = p.beta;
    var fwdMid = sqrt(p.fwd*K);
    var g1 = beta/fwdMid;
    return [g1,-g1*(1-beta)/fwdMid];
  }
  function theta(K){
    var beta = p.beta;
    return p.alpha/(p.sigma0*(1-beta))*(pow(p.fwd,1-beta)-pow(K,1-beta));
  }
  function normal(K){
    var alpha = p.alpha;
    var eps = alpha*alpha*p.T;
    var fwdMid = sqrt(p.fwd*K);
    var c = pow(fwdMid,p.beta);
    var g = gamma(K);
    var d = D(theta(K));
    var tmp = p.sigma0*c/alpha;
    var taylor = (1+eps*((2*g[1]-sq(g[0]))*sq(tmp)/24 + p.rho*g[0]*tmp/4 + (2-3*sq(p.rho))/24));
    if(abs(d) < 1e-8) return p.sigma0*pow(fwdMid,p.beta-1)*taylor;
    return alpha*(p.fwd-K)/d*taylor;
  }
  function logNormal(K){
    var alpha = p.alpha;
    var eps = alpha*alpha*p.T;
    var fwdMid = sqrt(p.fwd*K);
    var c = pow(fwdMid,p.beta);
    var g = gamma(K);
    var d = D(theta(K));
    var tmp = p.sigma0*c/alpha;
    var taylor = (1+eps*((2*g[1]-sq(g[0])+1/sq(fwdMid))*sq(tmp)/24 + p.rho*g[0]*tmp/4 + (2-3*sq(p.rho))/24));
    if(abs(d) < 1e-8) return p.sigma0*pow(fwdMid,p.beta-1)*taylor;
    return alpha*ln(p.fwd/K)/d*taylor;
  }
  var ret = isNormal?normal:logNormal;
  ret.getSigma0Normal = function(){
    var alpha = p.alpha;
    var eps = alpha*alpha*p.T;
    var fwdMid = p.fwd;
    var c = pow(fwdMid,p.beta);
    var g = gamma(p.fwd);
    var A = eps*(2*g[1]-sq(g[0]))*(c/alpha)*(c/alpha)/24;
    var B = eps*p.rho*g[0]*c/alpha/4;
    var C = 1 + eps*(2-3*sq(p.rho))/24;
    var D = -atmVol*pow(fwdMid,1-p.beta);
    return rootPolynome(A,B,C,D);
  }
  ret.transform = function(arr){
    p.sigma0 = exp(-arr[0]);
    p.alpha = exp(-arr[1]);
    p.beta = 1/(1+exp(-arr[2]));
    p.rho = 1-2/(1+exp(-arr[3]));
  }
  ret.invTransform = function(p){
    return [-ln(p.sigma0),-ln(p.alpha),-ln(1/p.beta-1),-ln(2/(1-p.rho)-1)];
  }
  ret.pts = [];
  ret.score = function(){
    var s = 0;
    ret.pts.forEach(function(pt){
      var tmp = ret(pt.k)-pt.v;
      s += tmp*tmp;
    });
    if(s < 0) s = -s;
    return s/ret.pts.length;
  }
  ret.points = function(arrT){
    return arrT.map(function(k){ return {k:k,v:ret(k)} });
  }
  ret.optimiseNum = function(first){
    var best = ret.score();
    var bestP = ret.invTransform(p);
    try {
      var guess;
      if(!first) guess = ret.invTransform({sigma0:rndUni(),alpha:rndUni(),rho:1-2*rndUni(),beta:rndUni()}); 
      else guess = ret.invTransform(p);
      var oracle = function(x){
        ret.transform(x);
        var s = ret.score();
        return s;
      }
      var res = numeric.uncmin(oracle,guess);
      if(res.f < best){
        best = res.f;
        bestP = res.solution;
      }
    }catch(ex){
      console.log("optimisation fail: ", ex);
    }
    ret.transform(bestP);
    return bestP;
  }
  ret.optimise = function(step){
    var nbTry = step;
    var first = true;
    while(nbTry--){
      ret.optimiseNum(first);
      first = false;
    }
  }
  return ret;
}

if (hasModule){
  module.exports.sabr = sabr;
  function test(){
    var p = {fwd:1,T:1,sigma0:0.1,alpha:0.1,beta:0.5,rho:-0.5};
    var s = sabr(p);
    s.pts.push({k:0.8,v:0.06});
    s.pts.push({k:0.9,v:0.055});
    s.pts.push({k:1.0,v:0.05});
    s.pts.push({k:1.1,v:0.048});
    s.pts.push({k:1.2,v:0.053});
    s.optimise(5);
    s.pts.forEach(function(pt){ console.log(pt.k,pt.v,s(pt.k)); });
  }
  if(require.main === module) return test();
}else globalScope.sabr = sabr;
}).call(this);