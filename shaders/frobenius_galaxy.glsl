// ======================================================================
//  Frobenius flow on the CM "galaxy"        — paste into Shadertoy's Image tab
// ======================================================================
//  Curve  y^2 = x^3 + 3x  over F_5.  Lattice  Lambda = Z[i].  Frobenius lifts
//  to multiplication by  alpha = -2 + i  on the torus  C/Lambda.  The points
//  are  E(F_125) = ker(alpha^3 - 1); the flow is  z -> alpha^phi * z, drawn in
//  the multiplicative ("galaxy") model  C^x / q^Z  via  w = exp(2 pi i z).
//
//  Shadertoy supplies iResolution, iTime and calls mainImage() for every pixel.
// ----------------------------------------------------------------------
#define N 130            // number of points = |alpha^3 - 1|^2
#define VIEW_GALAXY 1    // 0 = the flat torus instead

const float PI   = 3.14159265;
const vec2  ALPHA= vec2(-2.0, 1.0);
const float ARGA = 2.6779451;               // arg(alpha)
const float ABSA = 2.2360680;               // |alpha| = sqrt 5
const vec2  INVB = vec2(-3.0, -11.0)/130.0; // 1/(alpha^3 - 1) = conj(beta)/|beta|^2
const vec2  TAUP = vec2(0.4, 0.2);          // rebased basis tau' for the galaxy
const float QABS = 0.2846096;               // |q| = exp(-2 pi Im tau')

// complex numbers live in vec2:  z = x + i y  <->  vec2(x, y)
vec2 cmul(vec2 a, vec2 b){ return vec2(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x); }
vec2 apow(float p){ float r = pow(ABSA, p), t = ARGA*p; return r*vec2(cos(t), sin(t)); }
vec2 rnd(vec2 v){ return floor(v + 0.5); }  // nearest Gaussian integer (torus wrap)

vec3 hue(float ang){                        // angle (radians) -> colour wheel
  float h = fract(ang/(2.0*PI));
  vec3 p = abs(fract(h + vec3(0.0, 2.0/3.0, 1.0/3.0))*6.0 - 3.0);
  return mix(vec3(0.5), clamp(p - 1.0, 0.0, 1.0), 0.85);
}

// brightness of the flow at torus point z (basis 1, i) and phase phi in [0,1)
float field(vec2 z, float phi){
  float glow = 0.0;
  for(int k = 0; k < N; k++){
    vec2 g    = fract(float(k)*INVB);       // the k-th point,  k/beta  mod Z[i]
    vec2 prev = cmul(apow(phi), g);         // flowed head  alpha^phi * g
    prev     -= rnd(prev - z);              // pick the torus image nearest this pixel
    glow     += 1.5 * smoothstep(0.016, 0.0, length(z - prev));   // the dot
    for(int j = 1; j <= 12; j++){           // comet tail: short continuous segments
      vec2 cur = cmul(apow(phi - float(j)*0.0045), g);
      cur     -= rnd(cur - prev);           // keep the chain continuous
      vec2 pa  = z - prev, ba = cur - prev; // distance from z to segment [prev, cur]
      float d  = length(pa - ba*clamp(dot(pa, ba)/dot(ba, ba), 0.0, 1.0));
      float f  = 1.0 - float(j)/13.0;
      glow    += f*f * 0.7 * smoothstep(0.011, 0.0, d);
      prev = cur;
    }
  }
  return min(glow, 1.6);
}

void mainImage(out vec4 fragColor, in vec2 fragCoord){
  vec2 uv  = (fragCoord - 0.5*iResolution.xy)/iResolution.y;  // centre; y in [-.5,.5]
  float phi = fract(iTime*0.08);            // one Frobenius beat every 1/0.08 s
  vec3 col = vec3(0.055, 0.06, 0.067);
#if VIEW_GALAXY
  float R = length(uv), TH = atan(uv.y, uv.x), RAD = 0.46;
  float rn = R/RAD;
  if(rn <= 1.0 && rn >= QABS){              // inside the annulus
    float vp = -log(rn)/(2.0*PI*TAUP.y);    // radius -> vp    (log-polar = the galaxy map)
    float up = TH/(2.0*PI) - vp*TAUP.x;     // angle  -> up
    vec2  zp = vec2(up + vp*TAUP.x, vp*TAUP.y);      // z' = up*1 + vp*tau'
    col += hue(TH) * field(cmul(ALPHA, zp), phi);    // z = alpha*z'; colour by angle
  }
  col += vec3(0.10)*smoothstep(0.003, 0.0, abs(rn - 1.0));   // faint rim
#else
  vec2 z = uv/0.46;                         // flat torus, centred on 0
  if(max(abs(z.x), abs(z.y)) < 0.62) col += hue(atan(z.y, z.x)) * field(z, phi);
#endif
  fragColor = vec4(col, 1.0);
}
