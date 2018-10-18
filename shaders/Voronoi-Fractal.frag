/* Chaos City 2
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
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#version 130

uniform float iTime;
uniform vec2 iResolution;
uniform float iScale;
uniform float iNBeats;
uniform float iHighScale;

// Constants
const float pi = acos(-1.);
const vec3 c = vec3(1.,0.,-1.);

// Hash functions
float rand(vec2 x)
{
    return fract(sin(dot(x-1. ,vec2(12.9898,78.233)))*43758.5453);
}

float rand(vec3 x)
{
    return fract(sin(dot(x-1. ,vec3(12.9898,78.233,33.1818)))*43758.5453);
}

vec3 rand3(vec3 x)
{
    return vec3(rand(x.x*c.xx),rand(x.y*c.xx),rand(x.z*c.xx));
}

/* compute voronoi distance and closest point.
 * x: coordinate
 * return value: vec3(distance, coordinate of control point)
 */
vec3 vor(vec2 x)
{
    vec2 y = floor(x);
   	float ret = 1.;
    
    //find closest control point. ("In which cell am I?")
    vec2 pf=c.yy, p;
    float df=10., d;
    
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
        {
            p = y + vec2(float(i), float(j));
            p += vec2(rand(p), rand(p+1.));
            
            d = length(x-p);
            
            if(d < df)
            {
                df = d;
                pf = p;
            }
        }
    
    //compute voronoi distance: minimum distance to any edge
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
        {
            p = y + vec2(float(i), float(j));
            p += vec2(rand(p), rand(p+1.));
            
            vec2 o = p - pf;
            d = length(.5*o-dot(x-pf, o)/dot(o,o)*o);
            ret = min(ret, d);
        }
    
    return vec3(ret, pf);
}

// 2D value noise
float valuenoise(vec2 x)
{
    vec2 y = floor(x);
    x = fract(x);
    float r00 = -1.+2.*rand(y),
        r10 = -1.+2.*rand(y+c.xy),
        r01 = -1.+2.*rand(y+c.yx),
        r11 = -1.+2.*rand(y+c.xx);
    return mix(
        mix(r00, r10, smoothstep(0.,1.,x.x)),
        mix(r01, r11, smoothstep(0.,1.,x.x)),
        smoothstep(0.,1.,x.y)
    );
}

// 3D value noise
float valuenoise(vec3 x)
{
    vec3 y = floor(x);
    x = fract(x);
    float r000 = -1.+2.*rand(y),
        r100 = -1.+2.*rand(y+c.xyy),
        r010 = -1.+2.*rand(y+c.yxy),
        r001 = -1.+2.*rand(y+c.yyx),
        r110 = -1.+2.*rand(y+c.xxy),
        r011 = -1.+2.*rand(y+c.yxx),
        r101 = -1.+2.*rand(y+c.xyx),
        r111 = -1.+2.*rand(y+c.xxx);
    return 	mix(
        		mix(
            		mix(r000, r100, smoothstep(0.,1.,x.x)),
                    mix(r010, r110, smoothstep(0.,1.,x.x)),
                    smoothstep(0.,1.,x.y)
                ),
        		mix(
                    mix(r001, r101, smoothstep(0.,1.,x.x)),
                    mix(r011, r111, smoothstep(0.,1.,x.x)),
                    smoothstep(0.,1.,x.y)
                ),
        		smoothstep(0.,1.,x.z));
        
}

float mfvaluenoise(vec2 x, float fmin, float fmax, float phi)
{
    float sum = 0.;
    float a = 1.;
    
    for(float f = fmin; f<fmax; f = f*2.)
    {
        sum += a*valuenoise(f*x);
        a = a*phi;
    }
    
    return sum;
}

// add object to scene
vec2 add(vec2 sda, vec2 sdb)
{
    return mix(sda, sdb, step(sdb.x, sda.x));
}

// subtract object from scene
vec2 sub(vec2 sda, vec2 sdb)
{
    return mix(-sda, sdb, step(sda.x, sdb.x));
}

// Distance to line segment
float linesegment(vec3 x, vec3 p0, vec3 p1)
{
    vec3 d = p1-p0;
    float t = clamp(dot(x-p0,d)/dot(d,d),0.,1.);
    return length(x-mix(p0,p1,t));
}

// Stroke
float stroke(float sdf, float w)
{
    return abs(sdf)-w;
}

// extrusion
float zextrude(float z, float d2d, float h)
{
    vec2 d = abs(vec2(min(d2d, 0.),z))-h*c.yx;
    return min(max(d.x,d.y),0.)+length(max(d,0.));
}

/* compute voronoi distance and closest point.
 * x: coordinate
 * return value: vec3(distance, coordinate of control point)
 */
vec4 vor(vec3 x)
{
    vec3 y = floor(x);
   	float ret = 10.;
    
    //find closest control point. ("In which cell am I?")
    vec3 pf=c.yyy, p;
    float df=100., d;
    
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
            for(int k=-1; k<=1; k+=1)
            {
                p = y + vec3(float(i), float(j), float(k));
                p += rand3(p);

                d = length(x-p);
				
                if(d < df)
                {
                    df = d;
                    pf = p;
                }
            }
    
    //compute voronoi distance: minimum distance to any edge
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
            for(int k=-1; k<=1; k+=1)
            {
                p = y + vec3(float(i), float(j), float(k));
                p += rand3(p);

                vec3 o = p - pf;
                d = length(.5*o-dot(x-pf, o)/dot(o,o)*o);
                ret = min(ret, d);
            }
    return vec4(ret, pf);
}

mat3 rot(vec3 p)
{
    return mat3(c.xyyy, cos(p.x), sin(p.x), 0., -sin(p.x), cos(p.x))
        *mat3(cos(p.y), 0., -sin(p.y), c.yxy, sin(p.y), 0., cos(p.y))
        *mat3(cos(p.z), -sin(p.z), 0., sin(p.z), cos(p.z), c.yyyx);
}

float softmin( float a, float b, float k )
{
    float h = clamp( 0.5+0.5*(b-a)/k, 0.0, 1.0 );
    return mix( b, a, h ) - k*h*(1.0-h);
}

//BUILD A WORLD
//THEN BREAK IT
vec3 ind;
vec2 scene(vec3 x) // water
{
    x += iTime*c.yxy*1.e-1;
    
    vec3 dis = 12.*vec3((.1+.05*iScale)*valuenoise(x.xy-2.-2.e-1*iTime),(.1+.05*iScale)*valuenoise(x.yz-5.-2.e-1*iTime), (.1+.05*iScale)*valuenoise(x.zx-2.-2.e-1*iTime));
    vec4 v = vor(1.*x-dis);//, v2 = .1*vor(2.*x -2.*dis);
    ind = v.gba;
    float d = stroke(stroke(v.x, .1),.05);
    d = max(-d, -stroke(stroke(.5*v.x, .1),.01));
//     d *= v2.y;

    d = max(d, -length(x-c.yyx-iTime*c.yxy*1.e-1)+1.);

    d = min(d, length(x-mix(c.yxy, c.yyx, .8)-iTime*c.yxy*1.e-1)-.05-.3*iScale);
    d = softmin(d, length(max(abs(x-mix(c.yxy, c.yyx, .8)-.1*c.xyy-iTime*c.yxy*1.e-1)-(.05+.3*iScale)*c.xxx,0.)),max(.05-.55*clamp(iScale,0.,1.),.001));
// length(max(abs(p)-b,0.0))
    
//     float d = x.z - mfvaluenoise(x.xy-dis, 2., 40., .45+.2*clamp(3.*iScale, 0., 1.));
//     d = max(d, -.5*mfvaluenoise(x.xy-dis, 2., 10., .45+.2*clamp(3.*iScale,0., 1.)));
   
//     d = max(d, -mfvaluenoise(x.xy-dis, 50., 100., .45+.2*clamp(3.*iScale, 0.,1.)));
    //artificial guards for artifacts
    /*float dr = .465;
    vec3 y = mod(x, dr)-.5*dr;
    float guard = -length(max(abs(y)-vec3(.5*dr*c.xx, .6),0.));
    guard = abs(guard)+dr*.1;
    d = min(d, guard);
    */
    return vec2(d, 1.);
}

vec2 scene2(vec3 x) // tentacles
{
    x += iTime*c.yxy*1.e-2;
    
    vec2 dis = 12.*vec2((.1+.05*iScale)*valuenoise(x-2.-1.e-1*iTime),(.1+.05*iScale)*valuenoise(x.xy-5.-1.e-1*iTime));
    float d = x.z - mfvaluenoise(x.xy-dis, 6., 20., .45+.2*clamp(3.*iScale, 0., 1.));
    
    //artificial guards for artifacts
    float dr = .165;
    vec3 y = mod(x, dr)-.5*dr;
    float guard = -length(max(abs(y)-vec3(.5*dr*c.xx, .6),0.));
    guard = abs(guard)+dr*.1;
    d = min(d, guard);
    
    return vec2(d, 1.);
}

vec3 stdcolor(vec2 x)
{
	return 0.5 + 0.5*cos(iTime+x.xyx+vec3(0,2,4));
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

// Distance to line segment
float linesegment(vec2 x, vec2 p0, vec2 p1)
{
    vec2 d = p1-p0;
    float t = clamp(dot(x-p0,d)/dot(d,d),0.,1.);
    return length(x-mix(p0,p1,t));
}

// Distance to 210 logo
float logo(vec2 x, float r)
{
    return min(
        min(circle(x+r*c.zy, r), linesegment(x,r*c.yz, r*c.yx)),
        circlesegment(x+r*c.xy, r, -.5*pi, .5*pi)
    );
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
    u = c.yxx;
    t = c.yxy;
}

vec3 synthcol(float scale, float phase)
{
    vec3 c2 = vec3(207.,30.,102.)/255.,
        c3 = vec3(245., 194., 87.)/255.;
    mat3 r1 = rot((5.e-1*phase)*vec3(1.1,1.3,1.5));
    return 
        (
            1.1*mix
            (
                -(cross(c2, r1*c2)),
                -(r1*c2), 
                scale
            )
        );
}

vec3 color(float rev, float ln, float index, vec2 uv, vec3 x)
{
    vec3 col = c.yyy;
    if(index == 1.)
    {
   		vec3 c1 = stdcolor(x.xy+.5*rand(ind.xy+17.)+iNBeats), 
        	c2 = stdcolor(x.xy+x.yz+x.zx+.5*rand(ind.xy+12.)+iNBeats+11.+uv), 
            c3 = stdcolor(x.xy+x.yz+x.zx+.5*rand(ind.xy+15.)+iNBeats+23.+uv);
		col = .1*c1*vec3(1.,1.,1.) + .2*c1*vec3(1.,1.,1.)*ln + vec3(1.,1.,1.)*pow(rev,2.*(2.-1.5*clamp(iScale,0.,1.))) + 2.*c1*pow(rev, 8.)+3.*c1*pow(rev, 16.);
		col = clamp(.33*col, 0., 1.);
        //col = abs(col);
	}
    else if(index == 2.)
    {
        return stdcolor(x.xy+.5*rand(ind.xy+17.)+iNBeats);
    }
    return col;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord/iResolution.yy-.5;
    vec3 col = c.yyy;
    
    //if(iTime < 1000.) //scene 1
    {
        //use raymarching in this scene
    	vec3 ro, r, u, t, x, dir;
    	camerasetup(camera1, ro, r, u, t, uv, dir);
    	
        float d = 0.;// -(ro.z-1.)/dir.z;
    
        bool hit;
        vec2 s;
        raymarch(scene, x, ro, d, dir, s, 300, 1.e-4, hit);
        if(hit == false)
        {
            post(col, uv);
            fragColor = vec4(col, 1.);
            return;
        }
    
        vec3 n;
        calcnormal(scene, n, 1.e-3, x);
    
        vec3 l = x+2.*c.yyx, re = normalize(reflect(-l,n)), v = normalize(x-ro);
        float rev = abs(dot(re,v)), ln = abs(dot(l,n));
        
        col = color(rev, ln, s.y, uv, x);
        
        for(float i = .7; i >= .3; i -= .2)
        {
            //reflections
            dir = normalize(reflect(dir, n));
//             dir = normalize(refract(dir, n, i));
            d = 5.e-1;
            ro = x;
            raymarch(scene, x, ro, d, dir, s, 50, 5.e-4, hit);
            if(hit == false)
            {
                post(col, uv);
                fragColor = vec4(col, 1.);
                break;
            }
            calcnormal(scene, n, 1.e-3, x);
            l = x+2.*c.yyx;
            re = normalize(reflect(-l,n)); 
            v = normalize(x-ro);
            rev = abs(dot(re,v));
            ln = abs(dot(l,n));

            col = mix(col, color(rev, ln, s.y, uv, x), i);
        }
        
        //Portability
        //col = clamp((tanh(.7*col)), 0., 1.);
        //fog
        col = mix(col, c.yyy, tanh(2.e-1*(abs(x.y+x.z))));
        
    }
    
	post(col, uv);
    
    fragColor = vec4(col,1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
