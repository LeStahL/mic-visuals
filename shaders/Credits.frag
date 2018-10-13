/*
 * 210 Logo
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

/*************************************/
/** Constants and utility functions **/
/*************************************/

// Constants
const float pi = acos(-1.);
const vec3 c = vec3(1.,0.,-1.);

// Image rescaling function
float rescale(float value, vec2 old, vec2 new)
{
    return mix(new.x, new.y, (value-old.x)/(old.y-old.x));
}

/**********************/
/** Color generation **/
/**********************/

// Standard shadertoy color
vec3 stdcolor(vec2 x)
{
	return 0.5 + 0.5*cos(iTime+x.xyx+vec3(0,2,4));
}

/***********/
/** Noise **/
/***********/

// Hash function
float rand(vec2 x)
{
    return fract(sin(dot(x-1. ,vec2(12.9898,78.233)))*43758.5453);
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
        mix(r00, r10, x.x),
        mix(r01, r11, x.x),
        x.y
    );
}

// Multi-frequency value noise
float mfvaluenoise(vec2 x, float f0, float f1, float phi)
{
    float sum = 0.;
    float a = 1.2;
    
    for(float f = f0; f<f1; f = f*2.)
    {
        sum = a*valuenoise(f*x) + sum;
        a = a*phi;
    }
    
    return sum;
}

/****************/
/** 2D Objects **/
/****************/

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

/**********************/
/** Scene management **/
/**********************/

// Add objects to scene with proper antialiasing
vec4 add(vec4 sdf, vec4 sda)
{
    return vec4(
        min(sdf.x, sda.x), 
        mix(sda.gba, sdf.gba, smoothstep(-1.5/iResolution.y, 1.5/iResolution.y, sda.x))
    );
}

/******************/
/** 2D operators **/
/******************/

// Add width to object
float stroke(float sdf, float w)
{
    return abs(sdf)-w;
}

// Add thickness effect to object
vec4 thick(vec2 x, vec4 sdf, vec2 n)
{
    for(int i=1; i<6; ++i)
		sdf = add(vec4(stroke(sdf.x*n.x*n.y*2.*valuenoise((3.+4.*iScale)*x-2.-1.*iTime-1.2), .01), 3.e-3/abs(sdf.x+.2*valuenoise(x-2.-1.*iTime))*stdcolor(x+c.xx*.3*float(i))), sdf); 
    return sdf;
}

//distance to quadratic bezier spline with parameter t
float dist(vec2 p0,vec2 p1,vec2 p2,vec2 x,float t)
{
    t = clamp(t, 0., 1.);
    return length(x-pow(1.-t,2.)*p0-2.*(1.-t)*t*p1-t*t*p2);
}

//minimum distance to quadratic bezier spline
float dsp(vec2 p0, vec2 p1, vec2 p2, vec2 x)
{
    //coefficients for 0 = t^3 + a * t^2 + b * t + c
    vec2 E = x-p0, F = p2-2.*p1+p0, G = p1-p0;
    vec3 ai = vec3(3.*dot(G,F), 2.*dot(G,G)-dot(E,F), -dot(E,G))/dot(F,F);

	//discriminant and helpers
    float tau = ai.x/3., p = ai.y-tau*ai.x, q = - tau*(tau*tau+p)+ai.z, dis = q*q/4.+p*p*p/27.;
    
    //triple real root
    if(dis > 0.) 
    {
        vec2 ki = -.5*q*c.xx+sqrt(dis)*c.xz, ui = sign(ki)*pow(abs(ki), c.xx/3.);
        return dist(p0,p1,p2,x,ui.x+ui.y-tau);
    }
    
    //three distinct real roots
    float fac = sqrt(-4./3.*p), arg = acos(-.5*q*sqrt(-27./p/p/p))/3.;
    vec3 t = c.zxz*fac*cos(arg*c.xxx+c*pi/3.)-tau;
    return min(
        dist(p0,p1,p2,x, t.x),
        min(
            dist(p0,p1,p2,x,t.y),
            dist(p0,p1,p2,x,t.z)
        )
    );
}

//minimum distance to linear bezier spline
float dsg(vec2 p0, vec2 p1, vec2 x)
{
    vec2 d = p1-p0;
    float t = clamp(dot(x-p0,d)/dot(d,d),0.,1.);
    return length(x-mix(p0,p1,t));
}

// Draw Geometry
vec4 geometry(vec2 uv)
{
    //vec4 sdf = vec4(stroke(stroke(logo(x, .2), .06),.01), 2.5*stdcolor(x*1.7));
    float d = 1.;
    {
    const vec2 lin[88] = vec2[88](vec2(-4.00e-01,0.00e+00),vec2(-4.00e-01,6.67e-02),vec2(-4.00e-01,6.67e-02),vec2(-3.83e-01,6.67e-02),vec2(-4.00e-01,0.00e+00),vec2(-3.83e-01,0.00e+00),vec2(-3.67e-01,1.67e-02),vec2(-3.67e-01,5.00e-02),vec2(-3.49e-01,6.67e-02),vec2(-3.32e-01,6.67e-02),vec2(-3.32e-01,6.67e-02),vec2(-3.32e-01,1.67e-02),vec2(-2.81e-01,6.67e-02),vec2(-2.64e-01,0.00e+00),vec2(-2.64e-01,0.00e+00),vec2(-2.47e-01,6.67e-02),vec2(-2.30e-01,1.67e-02),vec2(-2.30e-01,3.33e-02),vec2(-2.30e-01,5.00e-02),vec2(-2.30e-01,5.00e-02),vec2(-1.96e-01,0.00e+00),vec2(-1.79e-01,0.00e+00),vec2(-1.79e-01,3.33e-02),vec2(-1.62e-01,3.33e-02),vec2(-1.45e-01,3.33e-02),vec2(-1.45e-01,1.67e-02),vec2(-1.12e-01,3.33e-02),vec2(-1.12e-01,0.00e+00),vec2(-6.08e-02,0.00e+00),vec2(-6.08e-02,3.33e-02),vec2(-4.33e-02,1.67e-02),vec2(-4.33e-02,6.67e-02),vec2(-9.17e-03,0.00e+00),vec2(7.50e-03,0.00e+00),vec2(7.50e-03,3.33e-02),vec2(2.42e-02,3.33e-02),vec2(7.58e-02,0.00e+00),vec2(7.58e-02,6.67e-02),vec2(7.58e-02,0.00e+00),vec2(9.25e-02,0.00e+00),vec2(7.58e-02,3.33e-02),vec2(9.25e-02,3.33e-02),vec2(1.27e-01,0.00e+00),vec2(1.27e-01,3.33e-02),vec2(2.12e-01,3.33e-02),vec2(2.12e-01,1.67e-02),vec2(2.45e-01,3.33e-02),vec2(2.45e-01,0.00e+00),vec2(2.96e-01,3.33e-02),vec2(2.96e-01,-1.67e-02),vec2(2.79e-01,-3.33e-02),vec2(2.63e-01,-3.33e-02),vec2(3.13e-01,0.00e+00),vec2(3.13e-01,6.67e-02),vec2(3.47e-01,0.00e+00),vec2(3.47e-01,1.67e-02),vec2(3.30e-01,3.33e-02),vec2(3.13e-01,3.33e-02),vec2(3.81e-01,1.67e-02),vec2(3.81e-01,6.67e-02),vec2(3.64e-01,3.33e-02),vec2(3.98e-01,3.33e-02),vec2(4.66e-01,1.67e-02),vec2(4.66e-01,6.67e-02),vec2(4.49e-01,3.33e-02),vec2(4.83e-01,3.33e-02),vec2(5.85e-01,3.33e-02),vec2(5.85e-01,1.67e-02),vec2(6.18e-01,3.33e-02),vec2(6.18e-01,-1.67e-02),vec2(6.02e-01,-3.33e-02),vec2(5.85e-01,-3.33e-02),vec2(6.87e-01,3.33e-02),vec2(6.87e-01,1.67e-02),vec2(7.20e-01,3.33e-02),vec2(7.20e-01,0.00e+00),vec2(7.72e-01,0.00e+00),vec2(7.72e-01,6.67e-02),vec2(7.72e-01,0.00e+00),vec2(7.88e-01,0.00e+00),vec2(7.72e-01,3.33e-02),vec2(7.88e-01,3.33e-02),vec2(8.22e-01,3.33e-02),vec2(8.22e-01,1.67e-02),vec2(8.56e-01,3.33e-02),vec2(8.56e-01,-1.67e-02),vec2(8.39e-01,-3.33e-02),vec2(8.22e-01,-3.33e-02)),
quad[150] = vec2[150](vec2(-3.83e-01,6.67e-02),vec2(-3.67e-01,6.67e-02),vec2(-3.67e-01,5.00e-02),vec2(-3.83e-01,0.00e+00),vec2(-3.67e-01,0.00e+00),vec2(-3.67e-01,1.67e-02),vec2(-3.32e-01,1.67e-02),vec2(-3.32e-01,0.00e+00),vec2(-3.49e-01,0.00e+00),vec2(-2.30e-01,1.67e-02),vec2(-2.30e-01,0.00e+00),vec2(-2.13e-01,0.00e+00),vec2(-1.79e-01,0.00e+00),vec2(-1.46e-01,8.33e-03),vec2(-1.79e-01,1.67e-02),vec2(-1.79e-01,1.67e-02),vec2(-2.12e-01,2.50e-02),vec2(-1.79e-01,3.33e-02),vec2(-1.45e-01,1.67e-02),vec2(-1.45e-01,0.00e+00),vec2(-1.28e-01,0.00e+00),vec2(-1.28e-01,0.00e+00),vec2(-1.12e-01,0.00e+00),vec2(-1.12e-01,1.67e-02),vec2(-7.75e-02,3.33e-02),vec2(-9.42e-02,3.33e-02),vec2(-9.42e-02,1.67e-02),vec2(-9.42e-02,1.67e-02),vec2(-9.42e-02,0.00e+00),vec2(-7.75e-02,0.00e+00),vec2(-7.75e-02,0.00e+00),vec2(-6.08e-02,0.00e+00),vec2(-6.08e-02,1.67e-02),vec2(-7.75e-02,3.33e-02),vec2(-6.08e-02,3.33e-02),vec2(-6.08e-02,1.67e-02),vec2(-4.33e-02,1.67e-02),vec2(-4.33e-02,0.00e+00),vec2(-2.67e-02,0.00e+00),vec2(7.50e-03,0.00e+00),vec2(4.08e-02,8.33e-03),vec2(7.50e-03,1.67e-02),vec2(7.50e-03,1.67e-02),vec2(-2.58e-02,2.50e-02),vec2(7.50e-03,3.33e-02),vec2(9.25e-02,0.00e+00),vec2(1.09e-01,0.00e+00),vec2(1.09e-01,1.67e-02),vec2(1.09e-01,1.67e-02),vec2(1.09e-01,3.33e-02),vec2(9.25e-02,3.33e-02),vec2(1.27e-01,1.67e-02),vec2(1.27e-01,3.33e-02),vec2(1.43e-01,3.33e-02),vec2(1.78e-01,0.00e+00),vec2(1.61e-01,0.00e+00),vec2(1.61e-01,1.67e-02),vec2(1.61e-01,1.67e-02),vec2(1.61e-01,3.33e-02),vec2(1.78e-01,3.33e-02),vec2(1.78e-01,3.33e-02),vec2(1.94e-01,3.33e-02),vec2(1.94e-01,1.67e-02),vec2(1.94e-01,1.67e-02),vec2(1.94e-01,0.00e+00),vec2(1.78e-01,0.00e+00),vec2(2.12e-01,1.67e-02),vec2(2.12e-01,0.00e+00),vec2(2.28e-01,0.00e+00),vec2(2.28e-01,0.00e+00),vec2(2.45e-01,0.00e+00),vec2(2.45e-01,1.67e-02),vec2(2.96e-01,1.67e-02),vec2(2.96e-01,3.33e-02),vec2(2.79e-01,3.33e-02),vec2(2.79e-01,3.33e-02),vec2(2.63e-01,3.33e-02),vec2(2.63e-01,1.67e-02),vec2(2.63e-01,1.67e-02),vec2(2.63e-01,0.00e+00),vec2(2.79e-01,0.00e+00),vec2(2.79e-01,0.00e+00),vec2(2.96e-01,0.00e+00),vec2(2.96e-01,1.67e-02),vec2(2.96e-01,-1.67e-02),vec2(2.96e-01,-3.33e-02),vec2(2.79e-01,-3.33e-02),vec2(3.47e-01,1.67e-02),vec2(3.47e-01,3.33e-02),vec2(3.30e-01,3.33e-02),vec2(3.81e-01,1.67e-02),vec2(3.81e-01,0.00e+00),vec2(3.98e-01,0.00e+00),vec2(4.66e-01,1.67e-02),vec2(4.66e-01,0.00e+00),vec2(4.83e-01,0.00e+00),vec2(5.17e-01,0.00e+00),vec2(5.00e-01,0.00e+00),vec2(5.00e-01,1.67e-02),vec2(5.00e-01,1.67e-02),vec2(5.00e-01,3.33e-02),vec2(5.17e-01,3.33e-02),vec2(5.17e-01,3.33e-02),vec2(5.33e-01,3.33e-02),vec2(5.33e-01,1.67e-02),vec2(5.33e-01,1.67e-02),vec2(5.33e-01,0.00e+00),vec2(5.17e-01,0.00e+00),vec2(5.85e-01,1.67e-02),vec2(5.85e-01,0.00e+00),vec2(6.02e-01,0.00e+00),vec2(6.02e-01,0.00e+00),vec2(6.18e-01,0.00e+00),vec2(6.18e-01,1.67e-02),vec2(6.18e-01,-1.67e-02),vec2(6.18e-01,-3.33e-02),vec2(6.02e-01,-3.33e-02),vec2(6.53e-01,0.00e+00),vec2(6.36e-01,0.00e+00),vec2(6.36e-01,1.67e-02),vec2(6.36e-01,1.67e-02),vec2(6.36e-01,3.33e-02),vec2(6.53e-01,3.33e-02),vec2(6.53e-01,3.33e-02),vec2(6.69e-01,3.33e-02),vec2(6.69e-01,1.67e-02),vec2(6.69e-01,1.67e-02),vec2(6.69e-01,0.00e+00),vec2(6.53e-01,0.00e+00),vec2(6.87e-01,1.67e-02),vec2(6.87e-01,0.00e+00),vec2(7.03e-01,0.00e+00),vec2(7.03e-01,0.00e+00),vec2(7.20e-01,0.00e+00),vec2(7.20e-01,1.67e-02),vec2(7.88e-01,0.00e+00),vec2(8.05e-01,0.00e+00),vec2(8.05e-01,1.67e-02),vec2(8.05e-01,1.67e-02),vec2(8.05e-01,3.33e-02),vec2(7.88e-01,3.33e-02),vec2(8.22e-01,1.67e-02),vec2(8.22e-01,0.00e+00),vec2(8.39e-01,0.00e+00),vec2(8.39e-01,0.00e+00),vec2(8.56e-01,0.00e+00),vec2(8.56e-01,1.67e-02),vec2(8.56e-01,-1.67e-02),vec2(8.56e-01,-3.33e-02),vec2(8.39e-01,-3.33e-02));
for(int i=0; i<44.0;++i) d=min(d,dsg(lin[2*i], lin[2*i+1], uv));
for(int i=0; i<50.0; ++i) d=min(d,dsp(quad[3*i], quad[3*i+1], quad[3*i+2], uv));
}
vec4 sdf = vec4(stroke(d, .006), 2.5*stdcolor(uv*1.7+iNBeats));
d = 1.;
{
const vec2 lin[48] = vec2[48](vec2(-2.00e-01,-2.00e-01),vec2(-2.00e-01,-6.67e-02),vec2(-2.00e-01,-6.67e-02),vec2(-1.33e-01,-2.00e-01),vec2(-1.33e-01,-2.00e-01),vec2(-1.33e-01,-6.67e-02),vec2(-9.83e-02,-2.00e-01),vec2(-9.83e-02,-6.67e-02),vec2(-9.83e-02,-6.67e-02),vec2(-6.50e-02,-6.67e-02),vec2(-9.83e-02,-1.33e-01),vec2(-6.50e-02,-1.33e-01),vec2(-3.17e-02,-1.67e-01),vec2(-3.17e-02,-2.00e-01),vec2(7.00e-02,-6.67e-02),vec2(3.33e-03,-1.67e-01),vec2(3.33e-03,-1.67e-01),vec2(7.00e-02,-1.67e-01),vec2(7.00e-02,-1.33e-01),vec2(7.00e-02,-2.00e-01),vec2(1.73e-01,-2.00e-01),vec2(2.40e-01,-6.67e-02),vec2(3.77e-01,-2.00e-01),vec2(3.77e-01,-6.67e-02),vec2(3.43e-01,-6.67e-02),vec2(4.10e-01,-6.67e-02),vec2(4.78e-01,-2.00e-01),vec2(5.12e-01,-2.00e-01),vec2(4.45e-01,-1.67e-01),vec2(5.12e-01,-1.67e-01),vec2(6.13e-01,-2.00e-01),vec2(6.13e-01,-1.33e-01),vec2(6.48e-01,-2.00e-01),vec2(6.48e-01,-1.33e-01),vec2(7.15e-01,-2.00e-01),vec2(7.15e-01,-1.67e-01),vec2(7.82e-01,-2.00e-01),vec2(7.82e-01,-1.67e-01),vec2(8.17e-01,-2.00e-01),vec2(8.83e-01,-2.00e-01),vec2(9.18e-01,-1.00e-01),vec2(9.52e-01,-6.67e-02),vec2(9.52e-01,-6.67e-02),vec2(9.52e-01,-2.00e-01),vec2(9.87e-01,-1.67e-01),vec2(9.87e-01,-1.00e-01),vec2(1.05e+00,-1.67e-01),vec2(1.05e+00,-1.00e-01)),
quad[57] = vec2[57](vec2(-6.50e-02,-6.67e-02),vec2(-3.17e-02,-6.67e-02),vec2(-3.17e-02,-1.00e-01),vec2(-3.17e-02,-1.00e-01),vec2(-3.17e-02,-1.33e-01),vec2(-6.50e-02,-1.33e-01),vec2(-6.50e-02,-1.33e-01),vec2(-3.17e-02,-1.33e-01),vec2(-3.17e-02,-1.67e-01),vec2(5.12e-01,-1.67e-01),vec2(5.12e-01,-1.33e-01),vec2(4.78e-01,-1.33e-01),vec2(4.78e-01,-1.33e-01),vec2(4.45e-01,-1.33e-01),vec2(4.45e-01,-1.67e-01),vec2(4.45e-01,-1.67e-01),vec2(4.45e-01,-2.00e-01),vec2(4.78e-01,-2.00e-01),vec2(5.80e-01,-1.33e-01),vec2(5.47e-01,-1.33e-01),vec2(5.47e-01,-1.67e-01),vec2(5.47e-01,-1.67e-01),vec2(5.47e-01,-2.00e-01),vec2(5.80e-01,-2.00e-01),vec2(5.80e-01,-2.00e-01),vec2(6.13e-01,-2.00e-01),vec2(6.13e-01,-1.67e-01),vec2(5.80e-01,-1.33e-01),vec2(6.13e-01,-1.33e-01),vec2(6.13e-01,-1.67e-01),vec2(6.48e-01,-1.67e-01),vec2(6.48e-01,-1.33e-01),vec2(6.82e-01,-1.33e-01),vec2(6.82e-01,-1.33e-01),vec2(7.15e-01,-1.33e-01),vec2(7.15e-01,-1.67e-01),vec2(7.15e-01,-1.67e-01),vec2(7.15e-01,-1.33e-01),vec2(7.48e-01,-1.33e-01),vec2(7.48e-01,-1.33e-01),vec2(7.82e-01,-1.33e-01),vec2(7.82e-01,-1.67e-01),vec2(8.17e-01,-6.67e-02),vec2(9.17e-01,-3.33e-02),vec2(8.17e-01,-2.00e-01),vec2(9.87e-01,-1.67e-01),vec2(9.87e-01,-2.00e-01),vec2(1.02e+00,-2.00e-01),vec2(1.02e+00,-2.00e-01),vec2(1.05e+00,-2.00e-01),vec2(1.05e+00,-1.67e-01),vec2(1.05e+00,-1.00e-01),vec2(1.05e+00,-6.67e-02),vec2(1.02e+00,-6.67e-02),vec2(1.02e+00,-6.67e-02),vec2(9.87e-01,-6.67e-02),vec2(9.87e-01,-1.00e-01));
for(int i=0; i<24.0;++i) d=min(d,dsg(lin[2*i], lin[2*i+1], uv));
for(int i=0; i<19.0; ++i) d=min(d,dsp(quad[3*i], quad[3*i+1], quad[3*i+2], uv));

}

    vec4 sda = vec4(stroke(d, .012), 2.5*stdcolor(uv*1.7+iNBeats+12.));
    //for(int i=0; i<10; ++i)
    //    sdf = add(sdf, vec4(stroke(circle(x-.5*vec2(valuenoise(x-2.-5.*iTime+2.*rand(float(i+3)*c.xx)), valuenoise(x-2.-5.*iTime+rand(float(i)*c.xx)))-.5*c.xy, .2+valuenoise(x-2.-5.*iTime+rand(float(i)*c.xx))),.01), 2.5*stdcolor(x+float(i)*.1)));
    return add(sdf,sda);
}

// Normal
const float dx = 1.e-4;
vec2 normal(vec2 x)
{
    float s = geometry(x).x;
    return normalize(vec2(geometry(x+dx*c.xy).x-s, geometry(x+dx*c.yx).x-s));
}

// Add Effects
vec4 scene(vec2 x)
{
    // Preprocessing
    x += .1*vec2(valuenoise(x-5.*iTime), valuenoise(x-2.-5.*iTime));
    
    // Normal scene
    vec4 sdf = geometry(x);
    vec2 n = normal(x);
    
    // Add potential effect
    sdf = thick(x, sdf, n);
    
    // Cover artefacts
    
    return sdf;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord/iResolution.yy-.5;
    vec4 s;
    //if(iTime < 20.)
    {
        s = scene(uv);
    }
    //else if(iTime < 40.)
    {
        
    }
    vec3 col = s.gba * smoothstep(1.5/iResolution.y, -1.5/iResolution.y, s.x);
    //small 210 logo
    col = mix(clamp(col,c.yyy,c.xxx), c.xxx, smoothstep(1.5/iResolution.y, -1.5/iResolution.y, stroke(logo(uv-vec2(-.45,.45),.02),.005)));
    
    col += vec3(0., 0.05, 0.1)*sin(uv.y*1050.+ 5.*iTime);

    fragColor = vec4(col,1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
