/* Puffin Pipes
 * From "Puffin Cure" by Team210 - 64k Demo at Under Construction 2018
 * Copyright (C) 2018  Alexander Kraus <nr4@z10.info>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#version 130

uniform float iTime;
uniform vec2 iResolution;
uniform float iScale;
uniform float iNBeats;
uniform float iHighScale;

// Global constants
const vec3 c = vec3(1.,0.,-1.);
const float pi = acos(-1.);

// Global variables
vec3 col = c.yyy;

// Hash function
float rand(vec2 x)
{
    return fract(sin(dot(x-1. ,vec2(12.9898,78.233)))*43758.5453);
}

/* Simplex noise -
Copyright (C) 2011 by Ashima Arts (Simplex noise)
Copyright (C) 2011-2016 by Stefan Gustavson (Classic noise and others)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
vec3 taylorInvSqrt(vec3 r) 
{     
    return 1.79284291400159-0.85373472095314*r; 
}

vec3 permute(vec3 x)
{
    return mod((x*34.+1.)*x, 289.);
}

float snoise(vec2 P) 
{     
    const vec2 C = vec2 (0.211324865405187134, 0.366025403784438597);  
    vec2 i = floor(P+dot(P, C.yy)) ; 
    vec2 x0 = P-i+dot(i, C.xx) ; 
    // Other  corners 
    vec2 i1 ; 
    i1.x = step ( x0.y , x0.x ) ;  //  1.0  i f  x0 . x > x0 . y ,  e l s e  0.0 
    i1.y = 1.0 - i1.x ; 
    // x1 = x0 − i1 + 1.0 ∗ C. xx ;  x2 = x0 − 1.0 + 2.0 ∗ C. xx ; 
    vec4 x12 = x0.xyxy + vec4 ( C.xx , C.xx * 2.0 - 1.0) ; 
    x12.xy -= i1 ; 
    //  Permutations 
    i = mod( i ,  289.0) ;  // Avoid  truncation  in  polynomial  evaluation 
    vec3 p = permute ( permute ( i.y + vec3 (0.0 , i1.y ,  1.0  ) ) + i.x + vec3 (0.0 , i1.x ,  1.0  ) ) ; 
    //  Circularly  symmetric  blending  kernel
    vec3 m = max(0.5 - vec3 ( dot ( x0 , x0 ) ,  dot ( x12.xy , x12.xy ) , dot ( x12.zw , x12.zw ) ) ,  0.0) ; 
    m = m * m ; 
    m = m * m ; 
    //  Gradients  from 41  points  on a  line ,  mapped onto a diamond 
    vec3 x = fract ( p * (1.0  /  41.0) ) * 2.0 - 1.0  ; 
    vec3 gy = abs ( x ) - 0.5  ; 
    vec3 ox = floor ( x + 0.5) ;  // round (x)  i s  a GLSL 1.30  feature 
    vec3 gx = x - ox ; //  Normalise  gradients  i m p l i c i t l y  by  s c a l i n g m 
    m *= taylorInvSqrt ( gx * gx + gy * gy ) ; // Compute  f i n a l  noise  value  at P 
    vec3 g ; 
    g.x = gx.x * x0.x + gy.x * x0.y ; 
    g.yz = gx.yz * x12.xz + gy.yz * x12.yw ; 
    //  Scale  output  to  span  range  [ − 1 ,1] 
    //  ( s c a l i n g  f a c t o r  determined by  experiments ) 
    return  -1.+2.*(130.0 * dot ( m , g ) ) ; 
}
/* End of Simplex Noise */

// Multi-frequency simplex noise
float mfsnoise(vec2 x, float f0, float f1, float phi)
{
    float sum = 0.;
    float a = 1.2;
    
    for(float f = f0; f<f1; f = f*2.)
    {
        sum = a*snoise(f*x) + sum;
        a = a*phi;
    }
    
    return sum;
}

// 3D rotational matrix
mat3 rot(vec3 p)
{
    return mat3(c.xyyy, cos(p.x), sin(p.x), 0., -sin(p.x), cos(p.x))
        *mat3(cos(p.y), 0., -sin(p.y), c.yxy, sin(p.y), 0., cos(p.y))
        *mat3(cos(p.z), -sin(p.z), 0., sin(p.z), cos(p.z), c.yyyx);
}

// add object to scene
vec2 add(vec2 sda, vec2 sdb)
{
    return mix(sda, sdb, step(sdb.x, sda.x));
}

// Distance to line segment
float lineseg(vec2 x, vec2 p1, vec2 p2)
{
    vec2 d = p2-p1;
    return length(x-mix(p1, p2, clamp(dot(x-p1, d)/dot(d,d),0.,1.)));
}

float lineseg(vec3 x, vec3 p1, vec3 p2)
{
    vec3 d = p2-p1;
    return length(x-mix(p1, p2, clamp(dot(x-p1, d)/dot(d,d),0.,1.)));
}

// Distance to circle
float circle(vec2 x, float r)
{
    return length(x)-r;
}

// Distance to circle segment
float circlesegment(vec2 x, float r, float p0, float p1)
{
    float p = atan(x.y, x.x);
    p = clamp(p, p0, p1);
    return length(x-r*vec2(cos(p), sin(p)));
}

// Distance to 210 logo
float logo(vec2 x, float r)
{
    return min(
        min(circle(x+r*c.zy, r), lineseg(x,r*c.yz, r*c.yx)),
        circlesegment(x+r*c.xy, r, -.5*pi, .5*pi)
    );
}

// Distance to stroke for any object
float stroke(float d, float w)
{
    return abs(d)-w;
}

vec2 scene(vec3 x)
{
    vec2 sdf = c.xy;
    
    x -= c.yyx*iTime;
    
    vec3 dv = 2.e-2*vec3(snoise(c.xx*(x.z-iTime)),
                      snoise(c.xx*(x.z-iTime)+13.),
                      snoise(c.xx*(x.z-iTime)-22.)
                      );
    x += dv;
    x *= rot(c.yyx*(iTime+iNBeats));
    
    float dz = .5;
    vec3 z = vec3(x.xy, mod(x.z, dz)-.5*dz);
    
    float r = length(x.xy),
        ddr = .1,
        dr = mod(r, ddr)-.5*ddr,
        p = atan(x.y,x.x),
        ddp = pi/(8.+round(24.*rand(iNBeats*c.xx+2.))),
        dp = mod(p, 2.*ddp)-ddp,
        r0 = .3;
    
    vec3 y = vec3(r0*cos(p-dp), r0*sin(p-dp), 0.);
    sdf = vec2(length(z-y)-.05, 1.+3.*round(rand(vec2(p-dp,x.z-z.z))));
    
    sdf = add(sdf, vec2(stroke(lineseg(z, vec3(r0*cos(p-dp-ddp), r0*sin(p-dp-ddp), 0.), vec3(r0*cos(p-dp+ddp), r0*sin(p-dp+ddp), 0.)), .03),1.));
	sdf = add(sdf, vec2(stroke(lineseg(vec3(z.xy,abs(z.z)), y, y+dz*c.yyx), .03), 1.+3.*round(rand(vec2(p-dp+17.)))));
    
    sdf.x -= length(dv)-.03*clamp(iScale,0.,1.);

    return sdf;
}

//performs raymarching
//scene: name of the scene function
//xc: 	 name of the coordinate variable
//ro:	 name of the ray origin variable
//d:	 name of the distance variable
//dir:	 name of the direction variable
//s:	 name of the scenestruct variable
//N:	 number of iterations used
//eps:	 exit criterion
//flag:  name of the flag to set if raymarching succeeded
#define raymarch(scene, xc, ro, d, dir, s, N, eps, flag) \
	flag = false;\
	for(int i=0; i<N; ++i)\
    {\
        xc = ro + d*dir;\
        s = scene(xc);\
        if(s.x < eps)\
        {\
            flag = true;\
            break;\
        }\
        d += s.x;\
    }

//computes normal with finite differences
//scene: name of the scene function
//n:	 name of the normal variable
//eps:	 precision of the computation
//xc:	 location of normal evaluation
#define calcnormal(scene, n, eps, xc) \
	{\
        float ss = scene(xc).x;\
        n = normalize(vec3(scene(xc+eps*c.xyy).xc-ss,\
                           scene(xc+eps*c.yxy).xc-ss,\
                           scene(xc+eps*c.yyx).xc-ss));\
    }

//camera setup
//camera: camera function with camera(out vec3 ro, out vec3 r, out vec3 u, out vec3 t)
//ro:	  name of the ray origin variable
//r:	  name of the right variable
//u:	  name of the up variable
//t:	  name of the target variable
//uv:	  fragment coordinate
//dir:	  name of the dir variable
#define camerasetup(camera, ro, r, u, t, uv, dir) \
	{\
        camera(ro, r, u, t);\
        t += uv.x*r+uv.y*u;\
        dir = normalize(t-ro);\
    }

//post processing: 210 logo and trendy display lines
//col: output color
//uv:  fragment coordinate
#define post(color, uv) \
	{\
    	col = mix(clamp(col,c.yyy,c.xxx), c.xxx, smoothstep(1.5/iResolution.y, -1.5/iResolution.y, stroke(logo(uv-vec2(-.45,.45),.02),.005)));\
    	col += vec3(0., 0.05, 0.1)*sin(uv.y*1050.+ 5.*iTime);\
	}
	
//camera for scene 1
void camera1(out vec3 ro, out vec3 r, out vec3 u, out vec3 t)
{
    ro = c.yyx;
    r = c.xyy;
    u = c.yxy;
    t = c.yyy;
}

// Standard shadertoy color
vec3 stdcolor(vec2 x)
{
	return 0.5 + 0.5*cos(iTime+x.xyx+vec3(0,2,4));
}

// Select material color
vec3 color(float rev, float ln, float index, vec2 uv, vec3 x)
{
    vec3 col = c.yyy;
    if(index == 1.)
    {
        x *= 1.e-2;
   		vec3 c1 = stdcolor(1.5e2*x.z+x.xy+.5*rand(17.*c.xx)+iNBeats), 
        	c2 = stdcolor(1.5e2*x.z+x.xy+x.yz+x.zx+.5*rand(12.*c.xx)+iNBeats+11.+uv), 
            c3 = stdcolor(1.5e2*x.z+x.xy+x.yz+x.zx+.5*rand(15.*c.xx)+iNBeats+23.+uv);
		col = .1*c1*vec3(1.,1.,1.) + .2*c1*vec3(1.,1.,1.)*ln + 1.5*vec3(1.,1.,1.)*pow(rev,2.*(2.)) + 2.*c1*pow(rev, 8.)+3.*c1*pow(rev, 16.);
        col = clamp(.23*col, 0., 1.);
	}
    else if(index == 4.)
    {
        x *= 1.e-2;
   		vec3 c1 = 1.*stdcolor(1.5e2*x.z+x.xy+.5*rand(47.*c.xx)+iNBeats+14.), 
        	c2 = 1.*stdcolor(1.5e2*x.z+x.xy+x.yz+x.zx+.5*rand(12.*c.xx)+iNBeats+21.+uv), 
            c3 = 1.*stdcolor(1.5e2*x.z+x.xy+x.yz+x.zx+.5*rand(15.*c.xx)+iNBeats+33.+uv);
		col = .1*c1*vec3(1.,1.,1.) + .2*c1*vec3(1.,1.,1.)*ln + 1.5*vec3(1.,1.,1.)*pow(rev,4.) + 2.*c1*pow(rev, 8.)+3.*c1*pow(rev, 16.);
        col = clamp(.23*col, 0., 1.);
	}
    return col;
}

float star(vec2 x, float r0)
{
    return 1.-smoothstep(.5*r0, r0, length(x));
}

// Background 1 (Moon)
vec3 background1(vec2 x)
{
    //Stars
    float dr = .03, scale;
    vec2 y = mod(x, dr)-.5*dr;
    float rs = rand(x-y)*.005,
        dx = -.5*(dr-rs)+(dr-2.*rs)*rand(x-y+1.),
        dy = -.5*(dr-rs)+(dr-2.*rs)*rand(x-y+2.);
    scale = star(y-vec2(dx,dy), rs);
    vec3 color = scale*clamp(8.*rand(x.xy+4.)*stdcolor(rand(x-y+3.)*x.xy), 0., 1.); 
    
    // Star nebula
    float f = mfsnoise(x.xy-6.93, 2.e-1, 1.e2, .55);
    color += mix(c.yyy, stdcolor(x), .5+.95*f);
    color += mix(c.yyy, 2.*stdcolor(x+4.), .5+.33*f);
    color += mix(c.yyy, stdcolor(x+8.), .5+.79*f);
    
    return clamp(color, 0., 1.);
}

// Background for the unc logo
vec3 background2(vec2 x)
{
    vec3 bg = c.yyy;
    float p = atan(x.y,x.x)/iTime,
        n = 5.,
        dmax = .3+.1*snoise(iTime*c.xx);
    for(float i = 0.; i<n; i+=1.)
    {
        float d = i/n*dmax;
        bg += background1((length(x)-.05+d-2.*iTime)*vec2(cos(p), sin(p))-.05*vec2(snoise(x.xy-iTime), snoise(x.xy+iTime)));
    }
    bg /= n;
    return bg;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord/iResolution.yy-.5;
    
	// Perform raymarching
    vec3 ro, r, u, t, x, dir;
    vec2 s;
    camerasetup(camera1, ro, r, u, t, uv, dir);
    	
    float d = .2/length(dir.xy);
    {
    	bool hit;
            
        raymarch(scene, x, ro, d, dir, s, 150, 1.e-4, hit);
        if(hit == false)
        {
        	// Draw Background here.
            col = background1(uv);
                
            post(col, uv);
            fragColor = vec4(col, 1.);
            return;
        }
            
        vec3 n;
        calcnormal(scene, n, 5.e-3, x);

        vec3 l = x+2.*c.yyx, re = normalize(reflect(-l,n)), v = normalize(x-ro);
        float rev = abs(dot(re,v)), ln = abs(dot(l,n));

        col = color(rev, ln, s.y, uv, x);
    }
    
    // Post-process
    post(col, uv);
    
    // Set the fragment color
    fragColor = vec4(col, 1.);    
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
