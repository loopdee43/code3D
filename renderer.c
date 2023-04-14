#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <gint/gint.h>
#include <gint/display.h>
#include <gint/display-cg.h>
#define SWAP(x,y) do { (x)=(x)^(y); (y)=(x)^(y); (x)=(x)^(y); } while(0)
typedef struct vec3d
{
	float x, y, z;
}vec3d;

typedef struct triangle
{
	vec3d p[3];
}triangle;

typedef struct mesh
{
	triangle tris[12];
}mesh;

typedef struct mat4x4
{
	float m[4][4];
}mat4x4;

int ScreenHeight(){
    return 224;
}

int ScreenWidth(){
    return 396;
}

void dpixel2(int x, int y, int color){
    int index = 396 * y + x;
    gint_vram[index] = color
}



// Horizontal line
void lcd_hline(uint8_t x1, uint8_t x2, uint8_t y, uint8_t color) {
	if(x1>=x2) SWAP(x1,x2);
	for(;x1<=x2;x1++) dpixel2(x1,y,color);
}

void fillTriangle(uint8_t x0, uint8_t y0,uint8_t x1, uint8_t y1, uint8_t x2, uint8_t y2, uint8_t color) {
 	uint8_t a, b, y, last;
  	// Sort coordinates by Y order (y2 >= y1 >= y0)
  	if (y0 > y1) { SWAP(y0, y1); SWAP(x0, x1); }
  	if (y1 > y2) { SWAP(y2, y1); SWAP(x2, x1); }
  	if (y0 > y1) { SWAP(y0, y1); SWAP(x0, x1); }
  
  	if(y0 == y2) { // All on same line case
    	a = b = x0;
    	if(x1 < a)      a = x1;
    	else if(x1 > b) b = x1;
    	if(x2 < a)      a = x2;
    	else if(x2 > b) b = x2;
        lcd_hline(a, b, y0, color);
        return;
    }

    int8_t
        dx01 = x1 - x0,
        dy01 = y1 - y0,
        dx02 = x2 - x0,
        dy02 = y2 - y0,
        dx12 = x2 - x1,
        dy12 = y2 - y1;
    int16_t sa = 0, sb = 0;

    // For upper part of triangle, find scanline crossings for segment
    // 0-1 and 0-2.  If y1=y2 (flat-bottomed triangle), the scanline y
    // is included here (and second loop will be skipped, avoiding a /
    // error there), otherwise scanline y1 is skipped here and handle
    // in the second loop...which also avoids a /0 error here if y0=y
    // (flat-topped triangle)
    if(y1 == y2) last = y1;   // Include y1 scanline
    else         last = y1-1; // Skip it

    for(y=y0; y<=last; y++) {
        a   = x0 + sa / dy01;
        b   = x0 + sb / dy02;
        sa += dx01;
        sb += dx02;
        // longhand a = x0 + (x1 - x0) * (y - y0) / (y1 - y0)
        //          b = x0 + (x2 - x0) * (y - y0) / (y2 - y0)
        lcd_hline(a, b, y);
    }

    // For lower part of triangle, find scanline crossings for segment
    // 0-2 and 1-2.  This loop is skipped if y1=y2
    sa = dx12 * (y - y1);
    sb = dx02 * (y - y0);
    for(; y<=y2; y++) {
        a   = x1 + sa / dy12;
        b   = x0 + sb / dy02;
        sa += dx12;
        sb += dx02;
        // longhand a = x1 + (x2 - x1) * (y - y1) / (y2 - y1)
        //          b = x0 + (x2 - x0) * (y - y0) / (y2 - y0)
        lcd_hline(a, b, y, color);
    }
}






void MultiplyMatrixVector(vec3d i, vec3d *o, mat4x4 m)
{
    o->x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + m.m[3][0];
	o->y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + m.m[3][1];
	o->z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + m.m[3][2];
	float w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + m.m[3][3];

	if (w != 0.0f)
	{
		o->x /= w; o->y /= w; o->z /= w;
	}
}

int main()
{
    mesh meshCube;
    mat4x4 matProj;
	vec3d vCamera;
    float fTheta;
    
    // SOUTH
    meshCube.tris[0].p[0] = (vec3d){0.0f, 0.0f, 0.0f};
    meshCube.tris[0].p[1] = (vec3d){0.0f, 1.0f, 0.0f};
    meshCube.tris[0].p[2] = (vec3d){1.0f, 1.0f, 0.0f};
    meshCube.tris[1].p[0] = (vec3d){0.0f, 0.0f, 0.0f};
    meshCube.tris[1].p[1] = (vec3d){1.0f, 1.0f, 0.0f};
    meshCube.tris[1].p[2] = (vec3d){1.0f, 0.0f, 0.0f};
    
    // EAST
    meshCube.tris[2].p[0] = (vec3d){1.0f, 0.0f, 0.0f};
    meshCube.tris[2].p[1] = (vec3d){1.0f, 1.0f, 0.0f};
    meshCube.tris[2].p[2] = (vec3d){1.0f, 1.0f, 1.0f};
    meshCube.tris[3].p[0] = (vec3d){1.0f, 0.0f, 0.0f};
    meshCube.tris[3].p[1] = (vec3d){1.0f, 1.0f, 1.0f};
    meshCube.tris[3].p[2] = (vec3d){1.0f, 0.0f, 1.0f};
    
    // NORTH
    meshCube.tris[4].p[0] = (vec3d){1.0f, 0.0f, 1.0f};
    meshCube.tris[4].p[1] = (vec3d){1.0f, 1.0f, 1.0f};
    meshCube.tris[4].p[2] = (vec3d){0.0f, 1.0f, 1.0f};
    meshCube.tris[5].p[0] = (vec3d){1.0f, 0.0f, 1.0f};
    meshCube.tris[5].p[1] = (vec3d){0.0f, 1.0f, 1.0f};
    meshCube.tris[5].p[2] = (vec3d){0.0f, 0.0f, 1.0f};
    
    // WEST
    meshCube.tris[6].p[0] = (vec3d){0.0f, 0.0f, 1.0f};
    meshCube.tris[6].p[1] = (vec3d){0.0f, 1.0f, 1.0f};
    meshCube.tris[6].p[2] = (vec3d){0.0f, 1.0f, 0.0f};
    meshCube.tris[7].p[0] = (vec3d){0.0f, 0.0f, 1.0f};
    meshCube.tris[7].p[1] = (vec3d){0.0f, 1.0f, 0.0f};
    meshCube.tris[7].p[2] = (vec3d){0.0f, 0.0f, 0.0f};
    
    // TOP
    meshCube.tris[8].p[0] = (vec3d){0.0f, 1.0f, 0.0f};
    meshCube.tris[8].p[1] = (vec3d){0.0f, 1.0f, 1.0f};
    meshCube.tris[8].p[2] = (vec3d){1.0f, 1.0f, 1.0f};
    meshCube.tris[9].p[0] = (vec3d){0.0f, 1.0f, 0.0f};
    meshCube.tris[9].p[1] = (vec3d){1.0f, 1.0f, 1.0f};
    meshCube.tris[9].p[2] = (vec3d){1.0f, 1.0f, 0.0f};
    
    // BOTTOM
    meshCube.tris[10].p[0] = (vec3d){1.0f, 0.0f, 1.0f};
    meshCube.tris[10].p[1] = (vec3d){0.0f, 0.0f, 1.0f};
    meshCube.tris[10].p[2] = (vec3d){0.0f, 0.0f, 0.0f};
    meshCube.tris[11].p[0] = (vec3d){1.0f, 0.0f, 1.0f};
    meshCube.tris[11].p[1] = (vec3d){0.0f, 0.0f, 0.0f};
    meshCube.tris[11].p[2] = (vec3d){0.0f, 0.0f, 1.0f};
    
    // Projection Matrix
    float fNear = 0.1f;
	float fFar = 1000.0f;
	float fFov = 90.0f;
	float fAspectRatio = (float)ScreenHeight() / (float)ScreenWidth();
	float fFovRad = 1.0f / tanf(fFov * 0.5f / 180.0f * 3.14159f);

	matProj.m[0][0] = fAspectRatio * fFovRad;
	matProj.m[1][1] = fFovRad;
	matProj.m[2][2] = fFar / (fFar - fNear);
	matProj.m[3][2] = (-fFar * fNear) / (fFar - fNear);
	matProj.m[2][3] = 1.0f;
	matProj.m[3][3] = 0.0f;
	
	int i = 0;
	float tp1 = clock();
	float tp2 = clock();
	
	while(i > -1){
	    tp2 = clock();
	    float fElapsedTime = tp2 - tp1;
	    mat4x4 matRotZ, matRotX;
	    fTheta += 5.0f
	    
	    matRotZ.m[0][0] = cosf(fTheta);
	    matRotZ.m[0][1] = sinf(fTheta);
		matRotZ.m[1][0] = -sinf(fTheta);
		matRotZ.m[1][1] = cosf(fTheta);
		matRotZ.m[2][2] = 1;
		matRotZ.m[3][3] = 1;

		// Rotation X
		matRotX.m[0][0] = 1;
		matRotX.m[1][1] = cosf(fTheta * 0.5f);
		matRotX.m[1][2] = sinf(fTheta * 0.5f);
		matRotX.m[2][1] = -sinf(fTheta * 0.5f);
		matRotX.m[2][2] = cosf(fTheta * 0.5f);
		matRotX.m[3][3] = 1;
		
		clearevents(); // prepping keyboard
		if (keydown(KEY_MENU)){
		    gint_osmenu();
		}
		dclear(); // clear screen
		for(int j; j < 11; j++){ // draw triangles
		    triangle triProjected, triTranslated, triRotatedZ, triRotatedZX;
		    
		    MultiplyMatrixVector(meshCube.tris[j].p[0], &triRotatedZ.p[0], matRotZ);
			MultiplyMatrixVector(meshCube.tris[j].p[1], &triRotatedZ.p[1], matRotZ);
			MultiplyMatrixVector(meshCube.tris[j].p[2], &triRotatedZ.p[2], matRotZ);

			// Rotate in X-Axis
			MultiplyMatrixVector(triRotatedZ.p[0], &triRotatedZX.p[0], matRotX);
			MultiplyMatrixVector(triRotatedZ.p[1], &triRotatedZX.p[1], matRotX);
			MultiplyMatrixVector(triRotatedZ.p[2], &triRotatedZX.p[2], matRotX);

			// Offset into the screen
			triTranslated = triRotatedZX;
			triTranslated.p[0].z = triRotatedZX.p[0].z + 3.0f;
			triTranslated.p[1].z = triRotatedZX.p[1].z + 3.0f;
			triTranslated.p[2].z = triRotatedZX.p[2].z + 3.0f;
			
			// finding surface normal
			vec3d normal, line1, line2;
			line1.x = triTranslated.p[1].x - triTranslated.p[0].x;
			line1.y = triTranslated.p[1].y - triTranslated.p[0].y;
			line1.z = triTranslated.p[1].z - triTranslated.p[0].z;

			line2.x = triTranslated.p[2].x - triTranslated.p[0].x;
			line2.y = triTranslated.p[2].y - triTranslated.p[0].y;
			line2.z = triTranslated.p[2].z - triTranslated.p[0].z;

			normal.x = line1.y * line2.z - line1.z * line2.y;
			normal.y = line1.z * line2.x - line1.x * line2.z;
			normal.z = line1.x * line2.y - line1.y * line2.x;

			// It's normally normal to normalise the normal
			float l = sqrtf(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
			normal.x /= l; normal.y /= l; normal.z /= l;
			
			if(normal.x * (triTranslated.p[0].x - vCamera.x) + 
			   normal.y * (triTranslated.p[0].y - vCamera.y) +
			   normal.z * (triTranslated.p[0].z - vCamera.z) < 0.0f)
		    {
				// Project triangles from 3D --> 2D
				MultiplyMatrixVector(triTranslated.p[0], &triProjected.p[0], matProj);
				MultiplyMatrixVector(triTranslated.p[1], &triProjected.p[1], matProj);
				MultiplyMatrixVector(triTranslated.p[2], &triProjected.p[2], matProj);

				// Scale into view
				triProjected.p[0].x += 1.0f; triProjected.p[0].y += 1.0f;
				triProjected.p[1].x += 1.0f; triProjected.p[1].y += 1.0f;
				triProjected.p[2].x += 1.0f; triProjected.p[2].y += 1.0f;
				triProjected.p[0].x *= 0.5f * (float)ScreenWidth();
				triProjected.p[0].y *= 0.5f * (float)ScreenHeight();
				triProjected.p[1].x *= 0.5f * (float)ScreenWidth();
				triProjected.p[1].y *= 0.5f * (float)ScreenHeight();
				triProjected.p[2].x *= 0.5f * (float)ScreenWidth();
				triProjected.p[2].y *= 0.5f * (float)ScreenHeight();

				// Rasterize triangle
				fillTriangle(triProjected.p[0].x, triProjected.p[0].y, triProjected.p[1].x, triProjected.p[1].y, triProjected.p[2].x, triProjected.p[2].y, 0);
			}
		}
		dupdate();
	}
	
}
