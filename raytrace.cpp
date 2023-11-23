using namespace std;

#include <stdio.h>
#include <iostream>
#include <glm/glm.hpp>
#include <array>
#include <cmath>
#include <vector>
#include <fstream>

// Output in P3 format, a text file containing:
// P3
// ncolumns nrows
// Max colour value (for us, and usually 255)
// r1 g1 b1 r2 g2 b2 .....
void save_imageP3(int Width, int Height, char* fname, unsigned char* pixels) {
	FILE *fp;
	const int maxVal=255;
	
    cout << fname;

	printf("Saving image %s: %d x %d\n", fname,Width,Height);
	fp = fopen(fname,"w");
	if (!fp) {
		printf("Unable to open file '%s'\n",fname);
		return;
	}
	fprintf(fp, "P3\n");
	fprintf(fp, "%d %d\n", Width, Height);
	fprintf(fp, "%d\n", maxVal);
	
	int k = 0 ;
	for(int j = 0; j < Height; j++) {
		
		for( int i = 0 ; i < Width; i++)
		{
			fprintf(fp," %d %d %d", pixels[k],pixels[k+1],pixels[k+2]) ;
			k = k + 3 ;
		}
		fprintf(fp,"\n") ;
	}
	fclose(fp);
}

string vecToString(glm::vec3 vec) {
    return "<" + to_string(vec.x) + ", " + to_string(vec.y) + ", " + to_string(vec.z) + ">\n";
}

struct Sphere {
    string name;

    glm::vec3 pos;
    glm::vec3 scale;
    glm::vec3 color;

    float ambient;
    float diffuse;
    float specular;
    float reflect;
    int n; 
};

struct Light {
    string name;

    glm::vec3 pos;
    glm::vec3 color;
};

struct Ray {
    glm::vec3 v;
    glm::vec3 S;
    int depth = 0;
    Sphere sphere;
};

class SceneInfo {
    public:
        float near;
        float left, right, top, bottom;
        int nCols, nRows;
        vector<Sphere> spheres;
        vector<Light> lights;
        glm::vec3 background;
        glm::vec3 ambient;
        string output;

        SceneInfo() = default;

        SceneInfo(string testFile) {
            ifstream file(testFile);
            string line;
            while(getline(file, line)) {
                vector<string> split = splitLine(line);
                parseLine(split);
            }
        }
    
    private:
        vector<string> splitLine(string line) {
            vector<string> split;
            int size = line.size();
            int prev = 0;
            for(int i=0; i<size; i++) {
                if(line[i] == ' ') {
                    split.push_back(line.substr(prev, i-prev));
                    prev = i+1;
                }
            }

            split.push_back(line.substr(prev, size-prev));

            return split;
        }

        void parseLine(vector<string> split) {
            string id = split[0];
            if(id == "NEAR") {
                this->near = stof(split[1]);

            } else if(id == "LEFT") {
                this->left = stof(split[1]);

            } else if(id == "RIGHT") {
                this->right = stof(split[1]);

            } else if(id == "BOTTOM") {
                this->bottom = stof(split[1]);

            } else if(id == "TOP") {
                this->top = stof(split[1]);

            } else if(id == "RES") {
                this->nCols = stof(split[1]);
                this->nRows = stof(split[2]);

            } else if(id == "SPHERE") {
                Sphere newSphere;
                newSphere.name = split[1];
                newSphere.pos = glm::vec3(
                    stof(split[2]), stof(split[3]), stof(split[4])
                );
                newSphere.scale = glm::vec3(
                    stof(split[5]), stof(split[6]), stof(split[7])
                );
                newSphere.color = glm::vec3(
                    stof(split[8]), 
                    stof(split[9]), 
                    stof(split[10])
                );
                newSphere.ambient = stof(split[11]);
                newSphere.diffuse = stof(split[12]);
                newSphere.specular = stof(split[13]);
                newSphere.reflect = stof(split[14]);
                newSphere.n = stoi(split[15]);

                this->spheres.push_back(newSphere);

            } else if(id == "LIGHT") {
                Light newLight;
                newLight.name = split[1];
                newLight.pos = glm::vec3(
                    stof(split[2]), stof(split[3]), stof(split[4])
                );
                newLight.color = glm::vec3(
                    stof(split[5]), 
                    stof(split[6]), 
                    stof(split[7])
                );

                this->lights.push_back(newLight);

            } else if(id == "BACK") {
                this->background = glm::vec3(
                    stof(split[1]), 
                    stof(split[2]), 
                    stof(split[3])
                );

            } else if(id == "AMBIENT") {
                this->ambient = glm::vec3(
                    stof(split[1]), 
                    stof(split[2]), 
                    stof(split[3])
                );

            } else if(id == "OUTPUT") {
                this->output = split[1];
            }
        }
};

const int nCol = 600;
const int nRow = 600;

const int MAX_DEPTH = 2;

const glm::vec3 eye = glm::vec3(0.0, 0.0, 0.0);
const glm::vec3 nCam = glm::vec3(0.0, 0.0, 1.0);
const glm::vec3 u = glm::vec3(1.0, 0.0, 0.0);
const glm::vec3 v = glm::vec3(0.0, 1.0, 0.0);

SceneInfo scene = SceneInfo();

Ray closestIntersection(Ray r) {

    Sphere closest;
    float t_min = FLT_MAX;
    for(int i=0; i<scene.spheres.size(); i++) {
        Sphere s = scene.spheres[i];

        // Transform ray
        Ray r_t;
        r_t.S = (r.S-s.pos)/s.scale;
        r_t.v = r.v/s.scale;

        float B = glm::dot(r_t.S, r_t.v);
        float B_square = pow(B, 2);
        float A = pow(glm::length(r_t.v), 2);
        float C = pow(glm::length(r_t.S), 2) - 1;

        // No or one intersection
        if( B_square - (A*C) <= 0) {
            continue;
        }

        float quadratic_1 = (-B/A);
        float quadratic_2 = sqrt(B_square - (A*C))/A;

        float t_plus  = quadratic_1 + quadratic_2;
        float t_minus = quadratic_1 - quadratic_2;

        float t = min(t_plus, t_minus);

        if(t < 0.00001) {
            continue;
        }

        if(t < t_min) {
            closest = scene.spheres[i];
            t_min = t;
        }
    }

    Ray intersection;
    if(t_min == FLT_MAX) {
        intersection.depth = MAX_DEPTH+1;
        return intersection;
    }

    // Intersection and normal
    intersection.S = glm::vec3(r.S + (t_min*r.v));
    intersection.v = glm::vec3( glm::normalize( 
        (intersection.S-closest.pos)/(closest.scale*closest.scale)
    ));

    intersection.sphere = closest;

    return intersection;
}

glm::vec4 raytrace(Ray r) {

    if(r.depth > MAX_DEPTH) {
        return glm::vec4(0.0,0.0,0.0,0.0);
    }

    Ray intersection = closestIntersection(r);
    if(intersection.depth == MAX_DEPTH+1) {
        return glm::vec4(scene.background, 0.0);
    }

    // Ambient product
    glm::vec3 ambient = intersection.sphere.color *
                        intersection.sphere.ambient *
                        scene.ambient;

    // Loop through lights
    glm::vec3 diffuse = glm::vec3(0.0, 0.0, 0.0);
    glm::vec3 specular = glm::vec3(0.0, 0.0, 0.0);
    Sphere sphere = intersection.sphere;

    for(int i=0; i<scene.lights.size(); i++) {
        Light light = scene.lights[i];
        glm::vec3 L = glm::vec3( glm::normalize(light.pos - intersection.S) );
        glm::vec3 N = intersection.v;

        // Diffuse
        float lightDotNormal = max( glm::dot(L, N), 0.0f );
        diffuse += sphere.color * 
                   lightDotNormal * 
                   light.color * 
                   sphere.diffuse;

        // Specular
        glm::vec3 R = glm::reflect(-L, N);
        glm::vec3 V = -1.0f * glm::normalize(intersection.S);
        float specularShiny = pow( max(glm::dot(R, V), 0.0f), sphere.n );

        if( glm::dot(L, N) < 0.0) {
            specular += glm::vec3(0.0, 0.0, 0.0);
        } else {
            specular += light.color *
                        sphere.specular * 
                        specularShiny;
        }
    }

    glm::vec4 clocal = glm::vec4( glm::vec3(ambient + diffuse + specular), 1.0);

    // Reflected ray
    Ray r_re;
    r_re.S = intersection.S;
    r_re.v = (-2 * ( glm::dot(intersection.v, r.v) ) * intersection.v) + r.v;
    r_re.depth = r.depth + 1;

    glm::vec4 c_re = raytrace(r_re);
    
    if(c_re[3] == 0.0) {
        return clocal;
    }

    return clocal + (sphere.reflect * c_re);
}

int main(int argc, char *argv[]) { 
    string testFile = argv[1];
    scene = SceneInfo(testFile);

    float W = (scene.right - scene.left)/2;
    float H = (scene.top - scene.bottom)/2;

    unsigned char *pixels;
    pixels = new unsigned char [3* 2*scene.nCols* 2*scene.nRows];
    int k = 0;

    for(int row=scene.nRows-1; row>=0; row--) {
        for(int col=0; col<scene.nCols; col++) {
            float u_c = (W*2*col/scene.nCols)-W;
            float v_r = (H*2*row/scene.nRows)-H;

            Ray r;
            r.v = glm::vec3((-scene.near*nCam) + (u_c*u) + (v_r*v));
            r.S = eye;
            r.depth = 1;

            glm::vec3 colors = raytrace(r);
            pixels[k] = min(colors.x*255.0f, 255.0f);
			pixels[k+1] = min(colors.y*255.0f, 255.0f);
			pixels[k+2] = min(colors.z*255.0f, 255.0f);

            k = k+3;
        }
    }

    save_imageP3(scene.nCols, scene.nRows, scene.output.data(), pixels);

}

