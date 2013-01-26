#pragma once

#include "ofMain.h"
#include <cmath>

typedef unsigned char ubyte;

struct GrassFire {

	ubyte* input;
	ubyte* output;
	int label;
	int m;
	int w, h;
	GrassFire(int w, int h){
		input = new ubyte[w*h];

		output = new ubyte[w*h];
		this->w = w;
		this->h = h;
		label = 0;
		m = 10;


	}

	void startFire(const ubyte* in, int w, int h){
		label = 1;
		memcpy(input, in, w*h);
		memset(output, 0, w*h);
		for(int y = 1 + m; y < h - 1 - m; y++){
			for(int x = 1 + m; x < w - 1 - m; x++){
				int i = y*w+x;
				if(input[i] == 1){
					fire(x, y);
					label++;
				}
			}
		}

	}
	void fire(int x, int y){
		if(x < 0 || x > w || y < 0 || y > h){
			return;
		}
		int i = y*w+x;
		if(input[i] == 0){
			return;
		}
		if(input[i] == 1){
			input[i] = 0;
			output[i] = label;

			fire(x, y - 1);
			fire(x + 1, y);
			fire(x, y + 1);
			fire(x - 1, y);
		}


	}

	static void blobToImage(ubyte* blob, ubyte* out, int w, int h){
		memset(out, 0, w*h*3);
		for(int i = 0; i < w*h; i++){

			int k = i*3;
			if(blob[i] != 0){
				int j = blob[i] % 2;
				out[k+j] = 255;
			}
		}
	}

	~GrassFire(){
		delete [] input;
		delete [] output;
	}



};

struct myPoint {
	int x, y;
	myPoint():x(0),y(0){

	}

	myPoint(int x, int y):x(x), y(y){

	}

	void rotate(float rads, const myPoint& orig){
		myPoint disp(x - orig.x, y - orig.y);
		x = disp.x * cos(rads) + disp.y * sin(rads);
		y = disp.y * cos(rads) - disp.x * sin(rads);
		x += orig.x;
		y += orig.y;
	}

	void set(int x, int y){
		this->x = x;
		this->y = y;
	}
};

struct BoundingBox {
	myPoint max;
	myPoint min;

	BoundingBox(){
		max.set(-10000, -10000);
		min.set(10000, 10000);
	}

	void set(myPoint p){
		if(p.x > max.x) max.x = p.x;
		if(p.y > max.y) max.y = p.y;
		if(p.x < min.x) min.x = p.x;
		if(p.y < min.y) min.y = p.y;
	}

	void set(int x, int y){
		set(myPoint(x, y));
	}

	int width(){
		return max.x - min.x;
	}

	int height(){
		return max.y - min.y;
	}

	myPoint getCenter(){
		return myPoint(width() / 2 + min.x, height() / 2 + min.y);
	}

	double getRatio(){
		return double(height()) / width();
	}
};

struct BLOB {

	myPoint centroid;
	BoundingBox box;
	int area;
	int label;
	unsigned int perim;
	BLOB(){
		area = 0;
		perim = 0;
	}

	double getCircularity(){
		return (perim) / (2.0 * sqrt(M_PI * area));
	}

	double getCompactness(){
		return double(area) / double(box.width() * box.height());
	}
};


class testApp : public ofBaseApp{

public:

	void setup();
	void update();
	void draw();

	void keyPressed(int key);
	void keyReleased(int key);
	void mouseMoved(int x, int y );
	void mouseDragged(int x, int y, int button);
	void mousePressed(int x, int y, int button);
	void mouseReleased(int x, int y, int button);
	void windowResized(int w, int h);
	void dragEvent(ofDragInfo dragInfo);
	void gotMessage(ofMessage msg);

	void pointProcess();
	void neighborhoodProcess();
	void morphology();
	void blobAnalysis();
	void transformations();



	ofVideoGrabber 		vidGrabber;
	ofImage				image;
	ubyte*		pointInput;
	ubyte*		pointOutput;
	ubyte*		binInput;
	ubyte*		binOutput;
	ubyte*		binDisplay;
	ubyte*		blackPixels;
	ubyte*		blendPixels;
	GrassFire*			gf;
	BLOB blob;
	ofTexture			videoTexture1;
	ofTexture			videoTexture2;

	ofTexture			geomTexture;
	ubyte*	geomInput;
	ubyte*	geomOutput;

	int 				w;
	int 				h;


};
