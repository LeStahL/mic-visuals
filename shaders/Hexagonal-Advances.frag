/*
 * Hexagonal Advances
 * 
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

const float pi = acos(-1.);
const vec3 c = vec3(1.,0.,-1.);
vec2 ind;

// hash function
float rand(vec2 a0)
{
    return fract(sin(dot(a0.xy ,vec2(12.9898,78.233)))*43758.5453);
}

// value noise
float vnoise(vec2 x)
{
    vec2 y = floor(x);
    x = fract(x);
    float r00 = -1.+2.*rand(y),
        r10 = -1.+2.*rand(y+c.xy),
        r01 = -1.+2.*rand(y+c.yx),
        r11 = -1.+2.*rand(y+c.xx);
    return mix(
        mix(r00, r10, x.x),
        mix(r01, r11, x.x),
        x.y
    );
}

float mfvnoise(vec2 x, float f0, float f1, float phi)
{
    float sum = 0.;
    float a = 1.;
    
    for(float f = f0; f<f1; f = f*2.)
    {
        sum = a*vnoise(f*x) + sum;
        a = a*phi;
    }
    
    return sum;
}

// compute distance to regular star
float dstar(vec2 x, float N, vec2 R)
{
    float d = pi/N,
        p0 = acos(x.x/length(x)),
        p = mod(p0, d),
        i = mod(round((p-p0)/d),2.);
    x = length(x)*vec2(cos(p),sin(p));
    vec2 a = mix(R,R.yx,i),
    	p1 = a.x*c.xy,
        ff = a.y*vec2(cos(d),sin(d))-p1;
   	ff = ff.yx*c.zx;
    return dot(x-p1,ff)/length(ff);
}

float zextrude(float z, float d2d, float h)
{
    vec2 d = abs(vec2(min(d2d, 0.),z))-h*c.yx;
    return min(max(d.x,d.y),0.)+length(max(d,0.));
}

float cr(vec2 x, float r, float w)
{
    return abs(length(x)-r)-w;
}

float cs(vec2 x, float r0, float w, float p0, float p1)
{
    float r = length(x), p = acos(x.x/r)*step(0.,x.y)-acos(x.x/r)*step(x.y,0.);
    p = clamp(p, p0, p1);
    vec2 y = r0*vec2(cos(p), sin(p));
    return length(x-y)-w;
}

float b(vec2 x, vec2 a, vec2 b, float w)
{
    vec2 d = b-a;
    return length(x-mix(a, b, clamp(dot(x-a, d)/dot(d,d), 0., 1.)))-w;
}

vec2 rot(vec2 x, float p)
{
    return mat2(cos(p), sin(p), -sin(p), cos(p))*x;
}

mat3 rot(vec3 p)
{
    vec3 cp = cos(p), sp = sin(p);
    mat3 m = mat3(cp.y*cp.x, cp.x*sp.z+cp.z*sp.x*sp.y, sp.x*sp.z-cp.x*cp.z*sp.y, 
           -cp.y*sp.z, cp.x*cp.z-sp.x*sp.y*sp.z, cp.z*sp.x+cp.x*sp.y*sp.z, 
           sp.y, -cp.y*sp.x, cp.x*cp.y);
    return m;
}

vec2 scene(vec3 x)
{
    
    float dr = .02, d, dh, h;
    vec3 y = vec3(mod(x.xy, dr)-.5*dr, x.z), index = x-y;
    
    vec2 z = (index.xy-c.yx)*(.5);
    float sd = min(cr(z-.125*c.xy, .125, .04), cs(z+.125*c.xy, .125, .04, -pi/2., pi/2.));
    sd = min(sd, b(z, -.125*c.yx, .125*c.yx, .04));
    
    dh = .05*mfvnoise(index.xy-.5*iTime, 8., 100., .45)+.1*clamp(iScale,0.,1.)*rand(index.xy);
    h = .1+dh+(.1+iScale*.1)*step(sd,0.);
    d = zextrude(y.z, dstar(rot(y.xy, 5.*iTime), 5.+5.*round(rand(3.*index.xy)), vec2(.2+.1*rand(2.*index.xy),.4+.1*rand(index.xy))*dr), h)-.001;
    
    ind = index.xy;
    
    if((sd < 0.) || (x.z < .5))
    {
    	float guard = -length(max(abs(y)-vec3(.5*dr*c.xx, .45),0.));
        guard = abs(guard)+dr*.1;
        d = min(d, guard);
    }
    
    vec2 sdf = vec2(d, 1.+step(sd, 0.)), sda = vec2(x.z, 1.);
    return mix(sdf, sda, step(sda.x, sdf.x));
}

const float dx = 1.e-3;
vec3 normal(vec3 x)
{
    float s = scene(x).x;
    return normalize(vec3(scene(x+dx*c.xyy).x-s, scene(x+dx*c.yxy).x-s, scene(x+dx*c.yyx).x-s));
}

vec3 background(vec2 x)
{
    return c.yyy;
}

vec3 synthcol(float scale, float phase)
{
    vec3 c2 = .5*vec3(rand(phase*c.xx), rand(phase*c.xx+1.), rand(phase*c.xx+2.))+.5;
    mat3 r1 = rot((5.e-1*phase)*vec3(1.1,1.3,1.5));
    float sc = rand(phase*c.xx);
    if(abs(sc) < .2)
        sc = sign(sc)*.2;
    return 
        (
            sc*1.1*mix
            (
                mix(-(cross(c2, r1*c2)),c.yyy, .5*scale),
                mix(c.yyy, -(r1*c2), .5*scale), 
                scale
            )
        );
}

// Stroke
float stroke(float sdf, float w)
{
    return abs(sdf)-w;
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

vec3 stdcolor(vec2 x)
{
	return 0.5 + 0.5*cos(iTime+x.xyx+vec3(0,2,4));
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord/iResolution.yy-.5, s;
    vec2 x0 = uv;
    uv.x += .1*vnoise(uv+iTime);
    vec3 or = rot(c.yyx*iTime)*(c.yyx-.2*c.yxy)+.5*c.yxy, ta = c.yxy, r = c.xyy, u = cross(r, normalize(ta-or)),
        rt = ta+uv.x*r+uv.y*u, rd = normalize(rt-or), x, col;
    
    //Initialize rays with intersection distance of ray and plane above scene
    float d = -(or.z-.45)/rd.z;
    
   	//Interval approximation root solving
   	float D = -(or.z+.1)/rd.z, dm;

    for(int i=0; i<300; ++i)
    {
        x = or + d * rd;
        s = scene(x);
        
        if(s.x < 1.e-4) break;
        if(i==299)
        {
            col = c.xxx*smoothstep(1.5/iResolution.y, -1.5/iResolution.y, stroke(logo(x0-vec2(-.45,.45),.02),.005));
            fragColor = vec4(col, 1.);
            return;
        }
        
        d += s.x;
    }

    vec3 n = normal(x), l = 2.*c.yyx, re = normalize(reflect(-l, n)), v = normalize(or-x);
    vec3 c1 = stdcolor(uv+2.5*ind.x+iNBeats+round(4.*x.z)), 
            c2 = stdcolor(uv+3.5*ind.y+iNBeats+round(4.*x.z)), 
            c3 = stdcolor(uv+3.5*ind.x+iNBeats+round(4.*x.z));
    float rev = clamp(dot(re,v),-1.,1.), ln = abs(dot(l,n));
    if(s.y == 1.)
    {
        col = .1*c1*vec3(1.,.3,.3) + .2*c1*vec3(1.,.3,.3)*ln + vec3(1.,1.,.1)*12.*tanh(x.z-.2)*pow(rev,2.*(2.-1.5*clamp(iScale,0.,1.))) + 2.*c1*pow(rev, 8.)+3.*c1*pow(rev, 16.);
        col = abs(col);
    }
    else if(s.y == 2.)
    {
        l = c.xxx;
        c1 += .1*mfvnoise(x.xz, 1., 100., .45);
        col = .1*c1*vec3(1.,.3,.3) + .2*c1*vec3(1.,.3,.3)*ln + vec3(1.,1.,.1)*pow(rev,2.*(2.-1.5*clamp(iScale,0.,1.))) + 2.*c1*pow(rev, 8.)+3.*c1*pow(rev, 16.);
        col = abs(col);
    }
    
    

    //portability
    col = 2.*clamp(.33*col, 0., 1.);
//     col = mix(col, c1*vec3(1.,.3,.3), tanh(1.e-2*length(x-ro)));
    //210 logo
    col = mix(clamp(col,c.yyy,c.xxx), c.xxx, smoothstep(1.5/iResolution.y, -1.5/iResolution.y, stroke(logo(x0-vec2(-.45,.45),.02),.005)));
    //trendy display lines
    col += vec3(0., 0.05, 0.1)*sin(x0.y*1050.+ 5.*iTime);
    
    fragColor = vec4(col,1.0);
}




void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
