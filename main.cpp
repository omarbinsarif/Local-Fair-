#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#define STB_IMAGE_IMPLEMENTATION
#pragma warning(disable:4996)
#include <unordered_map>
#include "shader.h"
#include "camera.h"
#include "basic_camera.h"
#include "pointLight.h"
#include "cube.h"
#include "sphere.h"
#include "stb_image.h"
#include <iostream>
#include <ctime>

using namespace std;

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void processInput(GLFWwindow* window);
void drawCube(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 model, float r, float g, float b);
void getCurrentTime(int& hours, int& minutes, int& seconds);
glm::mat4 transform(float tx, float ty, float tz, float sx, float sy, float sz);


/// floor 
void floor(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& floor);
void fence(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& floor);
void stall1call(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& stall1roof, Cube& stallstick, Cube& stallbase, Cube& stallbase2);
void stall2call(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& stall1roof, Cube& stallstick, Cube& stallbase, Cube& stallbase2);


//stall**********************
void stall1(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& stall1roof, Cube& stallstick, Cube& stallbase, Cube& stallbase2);
void stall2(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& stall1roof, Cube& stallstick, Cube& stallbase, Cube& stallbase2);
void drawCircle(GLfloat x, GLfloat y, GLfloat z, GLfloat radius, GLint numberOfSides);
void drawLine(GLfloat x1, GLfloat y1, GLfloat z1, GLfloat x2, GLfloat y2, GLfloat z2);
void triangle(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& stall1roof, Cube& stallstick, Cube& stallbase, Cube& stallbase2);
// Global variable to store the rotation angle
void drawCubes(glm::vec3 position, GLfloat size);
GLfloat angle = 0.0f;
void drawFlywheel(Shader& lightingShader, glm::vec3 translation, glm::vec3 scale, float rotationAngle, glm::vec3 color, unsigned int& cubeVAO, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& stall1roof);
// ************************DRAWING ROOM****************************************
void drawFlywheelstall(Shader& lightingShader, glm::vec3 translation, glm::vec3 scale, float rotationAngle, glm::vec3 color, unsigned int& cubeVAO, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& stall1roof);
void drawFlywheelcall(glm::vec3 translation, glm::vec3 scale, float rotationAngle, glm::vec3 color, glm::mat4 alTogether, glm::mat4 alTogether2, Shader& lightingShaderWithTexture, unsigned int& cubeVAO);

void canopi1(Shader& lightingShader, glm::vec3 color, glm::mat4 altogether);
void cotton(Shader& lightingShader, glm::vec3 color, glm::mat4 altogether);
void polls(Shader& lightingShader, glm::vec3 color, glm::mat4 altogether);
void drawCylinder( Shader& lightingShader, glm::vec3 color,glm::mat4 altogether);
void canopi(float rotz, Shader& lightingShader, glm::vec3 color, glm::mat4 altogether);
void canopistall(unsigned int& cubeVAO,Shader& lightingShader, glm::vec3 color, glm::mat4 alTogether, Cube& tex);
void canopistall2(unsigned int& cubeVAO, Shader& lightingShader, glm::vec3 color, glm::mat4 alTogether, Cube& tex);
void canopistall3(unsigned int& cubeVAO, Shader& lightingShader, glm::vec3 color, glm::mat4 alTogether, Cube& tex);
void drawCircle(GLfloat x, GLfloat y, GLfloat z, GLfloat radius, GLint numberOfSides);
void drawLine(GLfloat x1, GLfloat y1, GLfloat z1, GLfloat x2, GLfloat y2, GLfloat z2);
void drawFlywheel(glm::vec3 translation, glm::vec3 scale, float rotationAngle, glm::vec3 color, glm::mat4 alTogether, glm::mat4 alTogether2, Shader& lightingShaderWithTexture, unsigned int& cubeVAO);
void drawCircleflat(GLfloat x, GLfloat y, GLfloat z, GLfloat radius, GLint numberOfSides);
void drawride(glm::vec3 translation, glm::vec3 scale, float rotationAngle, glm::vec3 color, glm::mat4 alTogether, glm::mat4 alTogether2, Shader& lightingShaderWithTexture, unsigned int& cubeVAO);

unsigned int loadTexture(char const* path, GLenum textureWrappingModeS, GLenum textureWrappingModeT, GLenum textureFilteringModeMin, GLenum textureFilteringModeMax);
void shaderActivate(Shader& shader);

// settings
const unsigned int SCR_WIDTH = 1000;
const unsigned int SCR_HEIGHT = 800;

// modelling transform
float rotateAngle_X = 0.0;
float rotateAngle_Y = 0.0;
float rotateAngle_Z = 0.0;
float rotateAxis_X = 0.0;
float rotateAxis_Y = 0.0;
float rotateAxis_Z = 1.0;
float translate_X = 0.0;
float translate_Y = 0.0;
float translate_Z = 0.0;
float scale_X = 1.0;
float scale_Y = 1.0;
float scale_Z = 1.0;

// camera
Camera camera(glm::vec3(0.0f, 1.1f, 5.2f));
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
bool firstMouse = true;

float eyeX = 0.0, eyeY = 1.0, eyeZ = 3.0;
float lookAtX = 0.0, lookAtY = 0.0, lookAtZ = 0.0;
glm::vec3 V = glm::vec3(0.0f, 1.0f, 0.0f);
BasicCamera basic_camera(eyeX, eyeY, eyeZ, lookAtX, lookAtY, lookAtZ, V);


float width = 20.0;
float length = 20.0;

glm::vec3 corner1 = glm::vec3(-width / 2 - 0.5f, 0.0 - 20.80, -length / 2 - 0.80);
glm::vec3 corner2 = glm::vec3(width / 2 - 0.5, 0.0 - 20.80, -length / 2 - 0.80);
glm::vec3 corner3 = glm::vec3(width / 2 - 0.5, 0.0 - 20.80, length / 2 - 0.80);
glm::vec3 corner4 = glm::vec3(-width / 2 - 0.5, 0.0 - 20.80, length / 2 - 0.80);

// positions of the point lights
glm::vec3 pointLightPositions[] = {
   /* glm::vec3(0.17f, 0.4f, -1.75f),
    glm::vec3(0.0f,  1.5f,  0.0f),
    glm::vec3(0.0f,  1000.0f,  0.0f),
    glm::vec3(0.0f,  3.0f,  0.0f) 
    */

    corner1 + glm::vec3(0, 0.05, 0), // Add the height of the floor (0.05) to position the light above the floor
    corner2 + glm::vec3(0, 0.05, 0),
    corner3 + glm::vec3(0, 0.05, 0),
    corner4 + glm::vec3(0, 0.05, 0),
};






glm::vec3 point_light_positions[] = {
    glm::vec3(1.45f, 1.3f, 0.1f),
    glm::vec3(1.45f, 1.3f, -3.1f),
    glm::vec3(1.6f, 1.3f, -3.1f),
    glm::vec3(1.6f, 1.3f, 0.5f),
    glm::vec3(2.5f + 1.9f, 0.8f, -0.9f)

};

PointLight pointlight1(

    pointLightPositions[0].x, pointLightPositions[0].y, pointLightPositions[0].z,  // position
    0.7f, 0.7f, 0.7f,     // ambient
    0.7f, 0.7f, 0.7f,      // diffuse
    0.7f, 0.7f, 0.7f,        // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    1       // light number
);
PointLight pointlight2(

    pointLightPositions[1].x, pointLightPositions[1].y, pointLightPositions[1].z,  // position
    0.7f, 0.7f, 0.7f,     // ambient
    0.7f, 0.7f, 0.7f,      // diffuse
    0.7f, 0.7f, 0.7f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    2       // light number
);

PointLight pointlight3(

    pointLightPositions[2].x, pointLightPositions[2].y, pointLightPositions[2].z,  // position
    0.1f, 0.1f, 0.1f,     // ambient
    0.1f, 0.1f, 0.1f,      // diffuse
    0.1f, 0.1f, 0.1f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    3       // light number
);
PointLight pointlight4(

    pointLightPositions[3].x, pointLightPositions[3].y, pointLightPositions[3].z,  // position
    0.7f, 0.7f, 0.7f,     // ambient
    0.7f, 0.7f, 0.7f,      // diffuse
    0.7f, 0.7f, 0.7f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    4       // light number
);
// ******************************DRAWING_ROOM_LIGHT***********************************
PointLight light1(
    point_light_positions[0].x, point_light_positions[0].y, point_light_positions[0].z,  // position
    0.3f, 0.3f, 0.3f,     // ambient
    0.3f, 0.3f, 0.3f,      // diffuse
    0.3f, 0.3f, 0.3f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    5       // light number
);
PointLight light2(
    point_light_positions[1].x, point_light_positions[1].y, point_light_positions[1].z,  // position
    0.3f, 0.3f, 0.3f,     // ambient
    0.3f, 0.3f, 0.3f,      // diffuse
    0.3f, 0.3f, 0.3f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    6       // light number
);
PointLight light3(
    point_light_positions[2].x, point_light_positions[2].y, point_light_positions[2].z,  // position
    0.3f, 0.3f, 0.3f,     // ambient
    0.3f, 0.3f, 0.3f,      // diffuse
    0.3f, 0.3f, 0.3f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    7       // light number
);
PointLight light4(
    point_light_positions[3].x, point_light_positions[3].y, point_light_positions[3].z,  // position
    0.3f, 0.3f, 0.3f,     // ambient
    0.3f, 0.3f, 0.3f,      // diffuse
    0.3f, 0.3f, 0.3f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    8       // light number
);
PointLight light5(
    point_light_positions[4].x, point_light_positions[4].y, point_light_positions[4].z,  // position
    0.2f, 0.2f, 0.2f,     // ambient
    0.2f, 0.2f, 0.2f,      // diffuse
    0.2f, 0.2f, 0.2f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    9       // light number
);



// light settings
bool onOffPointToggle = true;
bool onOffSpotToggle = true;
bool onOffDirectToggle = true;
bool ambientToggle = true;
bool diffuseToggle = true;
bool specularToggle = true;

unsigned int ch_wood_tex,tex1;

bool isFlywheelRotating = false;
bool  startRotation1 = false;
float rtsp = 5.0f;
float flywheelRotationSpeed = 10.0f;
float flywheelRotation = 1;
bool startRotation = false;
//glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
//glm::mat4 view = camera.GetViewMatrix();
glm::mat4 projection;
glm::mat4 view;

string diffuseMapPath;
string specularMapPath;
GLuint texture0, texture1;
unsigned int tex2;

class Curve
{
public:
    vector<float> cntrlPoints;
    vector <float> coordinates;
    vector <float> normals;
    vector <int> indices;
    vector <float> vertices;
    vector<float> texCoords;
    const double pi = 3.14159265389;
    const int nt = 40;
    const int ntheta = 20;
    // Texture properties
    unsigned int diffuseMap;
    unsigned int specularMap;
    float shininess;
    Curve(vector<float>& tmp, unsigned int dMap, unsigned int sMap, float shiny)
        : diffuseMap(dMap), specularMap(sMap), shininess(shiny)
    {
        this->cntrlPoints = tmp;
        this->fishVAO = hollowBezier(cntrlPoints.data(), ((unsigned int)cntrlPoints.size() / 3) - 1);
        cout << cntrlPoints.size() << endl;
        cout << coordinates.size() << endl;
        cout << normals.size() << endl;
        cout << indices.size() << endl;
        cout << vertices.size() << endl;
    }
    ~Curve()
    {
        glDeleteVertexArrays(1, &fishVAO);
        glDeleteVertexArrays(1, &bezierVAO);
        glDeleteBuffers(1, &bezierVBO);
        glDeleteBuffers(1, &bezierEBO);
    }
    void draw(Shader& lightingShader, glm::mat4 model, glm::vec3 amb = glm::vec3(1.0f, 1.0f, 1.0f))
    {
        lightingShader.use();
        lightingShader.setMat4("model", model);
        lightingShader.setVec3("material.ambient", amb);
        lightingShader.setVec3("material.diffuse", amb);
        lightingShader.setVec3("material.specular", glm::vec3(0.5f, 0.5f, 0.5f));
        lightingShader.setFloat("material.shininess", 32.0f);

        // Set texture properties
        lightingShader.setInt("material.diffuseMap", 0);  // 0 corresponds to GL_TEXTURE0
        lightingShader.setInt("material.specularMap", 1); // 1 corresponds to GL_TEXTURE1

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, diffuseMap);

        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, specularMap);

        glBindVertexArray(fishVAO);
        glDrawElements(GL_TRIANGLES, (unsigned int)indices.size(), GL_UNSIGNED_INT, (void*)0);

        // unbind VAO
        glBindVertexArray(0);
    }
    void setTextureProperty(unsigned int dMap, unsigned int sMap, float shiny)
    {
        this->diffuseMap = dMap;
        this->specularMap = sMap;
        this->shininess = shiny;
    }


private:
    unsigned int fishVAO;
    unsigned int bezierVAO;
    unsigned int bezierVBO;
    unsigned int bezierEBO;


    unsigned int drawControlPoints()
    {
        unsigned int controlPointVAO;
        unsigned int controlPointVBO;

        glGenVertexArrays(1, &controlPointVAO);
        glGenBuffers(1, &controlPointVBO);

        glBindVertexArray(controlPointVAO);

        glBindBuffer(GL_ARRAY_BUFFER, controlPointVBO);
        glBufferData(GL_ARRAY_BUFFER, (unsigned int)cntrlPoints.size() * sizeof(float), cntrlPoints.data(), GL_STATIC_DRAW);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        return controlPointVAO;
    }

    long long nCr(int n, int r)
    {
        if (r > n / 2)
            r = n - r; // because C(n, r) == C(n, n - r)
        long long ans = 1;
        int i;

        for (i = 1; i <= r; i++)
        {
            ans *= n - r + i;
            ans /= i;
        }

        return ans;
    }
    void BezierCurve(double t, float xy[2], GLfloat ctrlpoints[], int L)
    {
        double y = 0;
        double x = 0;
        t = t > 1.0 ? 1.0 : t;
        for (int i = 0; i < L + 1; i++)
        {
            long long ncr = nCr(L, i);
            double oneMinusTpow = pow(1 - t, double(L - i));
            double tPow = pow(t, double(i));
            double coef = oneMinusTpow * tPow * ncr;
            x += coef * ctrlpoints[i * 3];
            y += coef * ctrlpoints[(i * 3) + 1];

        }
        xy[0] = float(x);
        xy[1] = float(y);
    }
    unsigned int hollowBezier(GLfloat ctrlpoints[], int L)
    {
        int i, j;
        float x, y, z, r;                // current coordinates
        float theta;
        float nx, ny, nz, lengthInv;    // vertex normal
        float u, v;                     // texture coordinates

        const float dtheta = 2 * pi / ntheta;        // angular step size

        float t = 0;
        float dt = 1.0 / nt;
        float xy[2];

        for (i = 0; i <= nt; ++i)              // step through y
        {
            BezierCurve(t, xy, ctrlpoints, L);
            r = xy[0];
            y = xy[1];
            theta = 0;
            t += dt;
            lengthInv = 1.0 / r;

            for (j = 0; j <= ntheta; ++j)
            {
                double cosa = cos(theta);
                double sina = sin(theta);
                z = r * cosa;
                x = r * sina;

                coordinates.push_back(x);
                coordinates.push_back(y);
                coordinates.push_back(z);

                // normalized vertex normal (nx, ny, nz)
                // center point of the circle (0,y,0)
                nx = (x - 0) * lengthInv;
                ny = (y - y) * lengthInv;
                nz = (z - 0) * lengthInv;

                normals.push_back(nx);
                normals.push_back(ny);
                normals.push_back(nz);

                // Calculate texture coordinates (s, t)
                u = static_cast<float>(j) / ntheta;
                v = static_cast<float>(i) / nt;
                texCoords.push_back(u);
                texCoords.push_back(v);

                theta += dtheta;
            }
        }
        // generate index list of triangles
        // k1--k1+1
        // |  / |
        // | /  |
        // k2--k2+1

        int k1, k2;
        for (int i = 0; i < nt; ++i)
        {
            k1 = i * (ntheta + 1);     // beginning of current stack
            k2 = k1 + ntheta + 1;      // beginning of next stack

            for (int j = 0; j < ntheta; ++j, ++k1, ++k2)
            {
                // k1 => k2 => k1+1
                indices.push_back(k1);
                indices.push_back(k2);
                indices.push_back(k1 + 1);

                // k1+1 => k2 => k2+1
                indices.push_back(k1 + 1);
                indices.push_back(k2);
                indices.push_back(k2 + 1);
            }
        }

        size_t count = coordinates.size();
        for (i = 0, j = 0; i < count; i += 3, j += 2)
        {
            //cout << count << ' ' << i + 2 << endl;
            vertices.push_back(coordinates[i]);
            vertices.push_back(coordinates[i + 1]);
            vertices.push_back(coordinates[i + 2]);

            if (i < normals.size())
                vertices.push_back(normals[i]);
            if (i + 1 < normals.size())
                vertices.push_back(normals[i + 1]);
            if (i + 2 < normals.size())
                vertices.push_back(normals[i + 2]);

            //// Add texture coordinates
            //if (j < texCoords.size())
            //    vertices.push_back(texCoords[j]);
            //if (j + 1 < texCoords.size())
            //    vertices.push_back(texCoords[j + 1]);
        }

        glGenVertexArrays(1, &bezierVAO);
        glBindVertexArray(bezierVAO);

        // create VBO to copy vertex data to VBO
        glGenBuffers(1, &bezierVBO);
        glBindBuffer(GL_ARRAY_BUFFER, bezierVBO);           // for vertex data
        glBufferData(GL_ARRAY_BUFFER,                   // target
            (unsigned int)vertices.size() * sizeof(float), // data size, # of bytes
            vertices.data(),   // ptr to vertex data
            GL_STATIC_DRAW);                   // usage

        // create EBO to copy index data
        glGenBuffers(1, &bezierEBO);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bezierEBO);   // for index data
        glBufferData(GL_ELEMENT_ARRAY_BUFFER,           // target
            (unsigned int)indices.size() * sizeof(unsigned int),             // data size, # of bytes
            indices.data(),               // ptr to index data
            GL_STATIC_DRAW);                   // usage

        // activate attrib arrays
        glEnableVertexAttribArray(0);
        glEnableVertexAttribArray(1);
        glEnableVertexAttribArray(2);
        // set attrib arrays with stride and offset
        int stride = 24;  // should be 24 bytes
        glVertexAttribPointer(0, 3, GL_FLOAT, false, stride, (void*)0);
        glVertexAttribPointer(1, 3, GL_FLOAT, false, stride, (void*)(sizeof(float) * 3));
        glVertexAttribPointer(2, 2, GL_FLOAT, false, stride, (void*)(sizeof(float) * 6)); // Add this line for texture coordinates
        //
                // unbind VAO, VBO and EBO
        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

        return bezierVAO;
    }

};

Curve* can;
Curve* can1;
Curve* can2;

vector<float>Bucket = {
0.001, 1.7054, 4.9697,
-0.550, 0.5803, 5.2077,
-1.00, -0.0164, 5.3340,
-1.500, -0.2268, 5.3785,

};


vector<float>Bucket1 = {
-1.50, 0.5803, 5.2077,
-1.00, -0.0164, 5.3340,
-0.900, -0.2268, 5.3785,
-0.00, -0.2268, 5.3785,

};

// timing
float deltaTime = 0.0f;    // time between current frame and last frame
float lastFrame = 0.0f;

int main()
{
    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // glfw window creation
    // --------------------
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "Local Fair", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);

    // tell GLFW to capture our mouse
    //glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    // configure global opengl state
    // -----------------------------
    glEnable(GL_DEPTH_TEST);

    // build and compile our shader zprogram
    // ------------------------------------
    Shader lightingShader("vertexShaderForPhongShading.vs", "fragmentShaderForPhongShading.fs");
    Shader lightingShaderWithTexture("vertexShaderForPhongShadingWithTexture.vs", "fragmentShaderForPhongShadingWithTexture.fs");
    Shader ourShader("vertexShader.vs", "fragmentShader.fs");

    // set up vertex data (and buffer(s)) and configure vertex attributes
    // ------------------------------------------------------------------

    float cube_vertices[] = {
        // positions      // normals
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        1.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        1.0f, 1.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        0.0f, 1.0f, 0.0f, 0.0f, 0.0f, -1.0f,

        1.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,
        1.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f,
        1.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f,
        1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f,

        0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f,
        1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f,
        1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f,
        0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f,

        0.0f, 0.0f, 1.0f, -1.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 1.0f, -1.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 0.0f, -1.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f, -1.0f, 0.0f, 0.0f,

        1.0f, 1.0f, 1.0f, 0.0f, 1.0f, 0.0f,
        1.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 1.0f, 1.0f, 0.0f, 1.0f, 0.0f,

        0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f,
        1.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f,
        1.0f, 0.0f, 1.0f, 0.0f, -1.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 0.0f, -1.0f, 0.0f
    };
    unsigned int cube_indices[] = {
        0, 3, 2,
        2, 1, 0,

        4, 5, 7,
        7, 6, 4,

        8, 9, 10,
        10, 11, 8,

        12, 13, 14,
        14, 15, 12,

        16, 17, 18,
        18, 19, 16,

        20, 21, 22,
        22, 23, 20
    };

    unsigned int cubeVAO, cubeVBO, cubeEBO;
    glGenVertexArrays(1, &cubeVAO);
    glGenBuffers(1, &cubeVBO);
    glGenBuffers(1, &cubeEBO);

    glBindVertexArray(cubeVAO);

    glBindBuffer(GL_ARRAY_BUFFER, cubeVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(cube_vertices), cube_vertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cubeEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(cube_indices), cube_indices, GL_STATIC_DRAW);


    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // vertex normal attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)12);
    glEnableVertexAttribArray(1);

    // second, configure the light's VAO (VBO stays the same; the vertices are the same for the light object which is also a 3D cube)
    unsigned int lightCubeVAO;
    glGenVertexArrays(1, &lightCubeVAO);
    glBindVertexArray(lightCubeVAO);

    glBindBuffer(GL_ARRAY_BUFFER, cubeVBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cubeEBO);
    // note that we update the lamp's position attribute's stride to reflect the updated buffer data
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);


    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);


    //stall1roof

    diffuseMapPath = "images/stall1roof.jpg";
    specularMapPath = "images/stall1roof.jpg";
    Cube stall1roof = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 5.0f, 5.0f, 6.0f, 6.0f);

    diffuseMapPath = "images/can.jpg";
    ch_wood_tex = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Curve Canopi(Bucket, ch_wood_tex, ch_wood_tex, 1.0f);
    can = &Canopi;

    texture0 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    texture1 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);

    diffuseMapPath = "images/stall2roof.jpg";
    tex1 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Curve Canopi1(Bucket, tex1, tex1, 1.0f);
    can1 = &Canopi1;

    texture0 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    texture1 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);

    diffuseMapPath = "images/cotton.png";
    tex2 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Curve Canopi2(Bucket1, tex2, tex2, 1.0f);
    can2 = &Canopi2;

    texture0 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    texture1 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);

   


    diffuseMapPath = "images/stallstick.jpg";
    specularMapPath = "images/stallstick.jpg";
    Cube stallstick = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 6.0f, 4.0f);


    diffuseMapPath = "images/stall2roof.jpg";
    specularMapPath = "images/stall2roof.jpg";
    Cube stallbase = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 5.0f, 5.0f, 6.0f, 6.0f);


    diffuseMapPath = "images/stallstick.jpg";
    specularMapPath = "images/stallstick.jpg";
    Cube stallbase2 = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 5.0f, 5.0f, 6.0f, 6.0f);


    diffuseMapPath = "images/curtain.png";
    specularMapPath = "images/curtain.png";
    Cube stall2back = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 5.0f, 5.0f, 6.0f, 6.0f);


    diffuseMapPath = "images/stallplaneroof.jpg";
    specularMapPath = "images/stallplaneroof.jpg";
    Cube stallplaneroof = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 5.0f, 5.0f, 6.0f, 6.0f);


    diffuseMapPath = "images/groun2.jpeg";
    specularMapPath = "images/groun2.jpeg";
    Cube ground = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 5.0f, 5.0f, 6.0f, 6.0f);


    diffuseMapPath = "images/brickwall.jpg";
    specularMapPath = "images/brickwall.jpg";
    Cube wall = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 5.0f, 5.0f, 6.0f, 6.0f);


    diffuseMapPath = "images/curtain.png";
    specularMapPath = "images/curtain.png";
    Cube curtain = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 5.0f, 5.0f, 6.0f, 6.0f);



   

    //pointlight2.turnOff();
    Sphere sphere;

    while (!glfwWindowShouldClose(window))
    {
        // per-frame time logic
               // --------------------
        float currentFrame = static_cast<float>(glfwGetTime());
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        // input
        // -----
        processInput(window);

        // render
        // ------
        glClearColor(0.298, 0.62, 0.961, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // be sure to activate shader when setting uniforms/drawing objects
        lightingShader.use();
        lightingShader.setVec3("viewPos", camera.Position);

        // activate shader
        lightingShader.use();

        projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
        lightingShader.setMat4("projection", projection);

        // camera/view transformation
        // glm::mat4 view = camera.GetViewMatrix();
        //glm::mat4 view = basic_camera.createViewMatrix();
        view = camera.GetViewMatrix();
        lightingShader.setMat4("view", view);



        // Modelling Transformation
        glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
        glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model,rotateMatrix;



        translateMatrix = glm::translate(identityMatrix, glm::vec3(translate_X, translate_Y, translate_Z));
        rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        scaleMatrix = glm::scale(identityMatrix, glm::vec3(scale_X, scale_Y, scale_Z));
        model = translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;
        lightingShaderWithTexture.setMat4("model", model);

        floor(cubeVAO, lightingShader, model, lightingShaderWithTexture, ground);
        lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1,1,1));

        translateMatrix = glm::translate(identityMatrix, glm::vec3(translate_X, translate_Y, translate_Z));
        rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        scaleMatrix = glm::scale(identityMatrix, glm::vec3(scale_X, scale_Y, scale_Z));
        model = translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;
        lightingShaderWithTexture.setMat4("model", model);

        fence(cubeVAO, lightingShader, model, lightingShaderWithTexture, stallstick);


       

        translateMatrix = glm::translate(identityMatrix, glm::vec3(translate_X, translate_Y, translate_Z));
        rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        scaleMatrix = glm::scale(identityMatrix, glm::vec3(scale_X, scale_Y, scale_Z));
        model = translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix ;

        lightingShaderWithTexture.setMat4("model", model);
        stall1call(cubeVAO, lightingShader, model, lightingShaderWithTexture, stall1roof, stallstick, stallbase, stallbase2);


        translateMatrix = glm::translate(identityMatrix, glm::vec3(translate_X, translate_Y, translate_Z));
        rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        scaleMatrix = glm::scale(identityMatrix, glm::vec3(scale_X, scale_Y, scale_Z));
        model = translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;

        lightingShaderWithTexture.setMat4("model", model);
        stall2call(cubeVAO, lightingShader, model, lightingShaderWithTexture, stall1roof, stallstick, stallbase, stallbase2);


        //translateMatrix = glm::translate(identityMatrix, glm::vec3(-1.5, 1, translate_Z));
        //rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        //rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        //rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        //scaleMatrix = glm::scale(identityMatrix, glm::vec3(.5, .5, .5));
        //model = translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;
       
        //lightingShaderWithTexture.setMat4("model", model);
        //glm::vec3 color = glm::vec3(7.0, 3.0, 1.0);
        //lightingShaderWithTexture.setVec3("curveColor", color);

        //
        //bucket->draw(lightingShaderWithTexture, model); // Set color to red
        //
        //lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1.0, 1.0, 1.0));


        translateMatrix = glm::translate(identityMatrix, glm::vec3(translate_X,translate_Y, translate_Z));
        rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        scaleMatrix = glm::scale(identityMatrix, glm::vec3(scale_X,scale_Y,scale_Z));
        model = translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;
        lightingShaderWithTexture.setMat4("model", model);

        canopistall(cubeVAO,lightingShaderWithTexture, glm::vec3(1.0f, 0.0f, 0.0f),model,curtain);

        lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1, 1, 1));


        translateMatrix = glm::translate(identityMatrix, glm::vec3(translate_X, translate_Y, translate_Z));
        rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        scaleMatrix = glm::scale(identityMatrix, glm::vec3(scale_X, scale_Y, scale_Z));
        model = translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;
        lightingShaderWithTexture.setMat4("model", model);

        //polls(lightingShaderWithTexture, glm::vec3(1.0f, 0.0f, 0.0f), model);

        lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1, 1, 1));


        translateMatrix = glm::translate(identityMatrix, glm::vec3(translate_X, translate_Y, translate_Z));
        rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        scaleMatrix = glm::scale(identityMatrix, glm::vec3(scale_X, scale_Y, scale_Z));
        model = translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;
        lightingShaderWithTexture.setMat4("model", model);

        canopistall3(cubeVAO, lightingShaderWithTexture, glm::vec3(1.0f, 0.0f, 0.0f), model, curtain);

        lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1, 1, 1));


        translateMatrix = glm::translate(identityMatrix, glm::vec3(translate_X, translate_Y, translate_Z));
        rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        scaleMatrix = glm::scale(identityMatrix, glm::vec3(scale_X, scale_Y, scale_Z));
        model = translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;
        lightingShaderWithTexture.setMat4("model", model);

        canopistall2(cubeVAO, lightingShaderWithTexture, glm::vec3(1.0f, 0.0f, 0.0f), model, curtain);

        lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1, 1, 1));


        translateMatrix = glm::translate(identityMatrix, glm::vec3(translate_X, translate_Y, translate_Z));
        rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        scaleMatrix = glm::scale(identityMatrix, glm::vec3(scale_X, scale_Y, scale_Z));
        model = translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;
        lightingShaderWithTexture.setMat4("model", model);

        cotton(lightingShaderWithTexture, glm::vec3(1.0f, 1.0f, 1.0f), model);

        lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1, 1, 1));

 /*       translateMatrix = glm::translate(identityMatrix, glm::vec3(translate_X, translate_Y, translate_Z));
        rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        scaleMatrix = glm::scale(identityMatrix, glm::vec3(1, .5, scale_Z));
        model = translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;
        lightingShaderWithTexture.setMat4("model", model);
        drawFlywheel(glm::vec3(0,0,1), model, lightingShaderWithTexture);*/



        if (isFlywheelRotating)
        {
            glm::vec3 flywheelPosition(translate_X, translate_Y, translate_Z);
            glm::vec3 flywheelScale(scale_X, scale_Y, scale_Z);

            glm::vec3 flywheelColor(0.0f, 0.0f, 10.0f);
            glm::mat4 model1, translateM;
            translateM = glm::mat4(1.0);
            translateM = glm::translate(translateM, glm::vec3(flywheelPosition));
            rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
            rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
            rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
           // glm::mat4 rotateFlywheelMatrix = glm::rotate(identityMatrix, glm::radians(flywheelRotation), glm::vec3(0.0f, 0.0f, 1.0f));
            model = rotateXMatrix * rotateYMatrix * rotateZMatrix ;
            model1 = translateM * rotateXMatrix * rotateYMatrix * rotateZMatrix;
            lightingShaderWithTexture.setMat4("model", model);


            drawFlywheel(flywheelPosition, flywheelScale, flywheelRotation, flywheelColor, model, model1, lightingShaderWithTexture, cubeVAO);


            flywheelRotation += 1;
        }

        else
        {
            //glm::vec3 flywheelPosition(-4.0f, 2.0f, -4.0f); // Position
            glm::vec3 flywheelPosition(translate_X,translate_Y,translate_Z);
            glm::vec3 flywheelScale(scale_X,scale_Y,scale_Z);     // Scale

            glm::vec3 flywheelColor(0.0f, 0.0f, 10.0f);

            glm::mat4 model1, translateM;
            translateM = glm::mat4(1.0);
            translateM = glm::translate(translateM, glm::vec3(flywheelPosition));
            rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
            rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
            rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
            model1 = translateM * rotateXMatrix * rotateYMatrix * rotateZMatrix;
            lightingShaderWithTexture.setMat4("model", model);


            drawFlywheel(flywheelPosition, flywheelScale, flywheelRotation, flywheelColor, model, model1, lightingShaderWithTexture, cubeVAO);



        }

        


      // glm::mat4 model1 = glm::mat4(1.0f);
     /*   glm::vec3 flywheelPosition(translate_X, translate_Y, translate_Z);
        glm::vec3 flywheelScale(scale_X,scale_Y,scale_Z);
        glm::vec3 flywheelColor(0,0,1);

        glm::mat4 model1, translateM;
        translateM = glm::mat4(1.0);
        translateM = glm::translate(translateM, glm::vec3(translate_X,translate_Y,translate_Z));
        rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        scaleMatrix = glm::scale(identityMatrix, glm::vec3(scale_X, scale_Y, scale_Z));
        model = translateM  * scaleMatrix  * rotateXMatrix * rotateYMatrix * rotateZMatrix;
        model1 = translateM *scaleMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix;
        lightingShaderWithTexture.setMat4("model", model);


        drawFlywheelcall(flywheelPosition, flywheelScale, flywheelRotation, flywheelColor, model, model1, lightingShaderWithTexture, cubeVAO);

*/
        glm::vec3 flywheelPosition(translate_X, translate_Y, translate_Z);;
        glm::vec3 flywheelScale(scale_X,scale_Y,scale_Z);     

        glm::vec3 flywheelColor(0.0f, 2.0f, 1.0f);

        glm::mat4 model1, translateM;
        translateM = glm::mat4(1.0);
        translateM = glm::translate(translateM, glm::vec3(flywheelPosition));
        rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        model =rotateXMatrix * rotateYMatrix * rotateZMatrix ;
        model1 = translateM * rotateXMatrix * rotateYMatrix * rotateZMatrix;
        lightingShaderWithTexture.setMat4("model", model);


        //drawtrain(flywheelPosition, flywheelScale, flywheelRotation, flywheelColor, model, model1, lightingShaderWithTexture, cubeVAO);

        drawride(flywheelPosition, flywheelScale, flywheelRotation, flywheelColor, model, model1, lightingShaderWithTexture, cubeVAO);
        lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1, 1, 1));
        

        lightingShader.use();
        // point light 1
        pointlight1.setUpPointLight(lightingShaderWithTexture);
        // point light 2
        pointlight2.setUpPointLight(lightingShaderWithTexture);
        // point light 3
        pointlight3.setUpPointLight(lightingShaderWithTexture);
        // point light 4
        pointlight4.setUpPointLight(lightingShaderWithTexture);

       
        light1.setUpPointLight(lightingShaderWithTexture);
        light2.setUpPointLight(lightingShaderWithTexture);
        light3.setUpPointLight(lightingShaderWithTexture);
        light4.setUpPointLight(lightingShaderWithTexture);
        light5.setUpPointLight(lightingShaderWithTexture);
      

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // optional: de-allocate all resources once they've outlived their purpose:
    // ------------------------------------------------------------------------
    glDeleteVertexArrays(1, &cubeVAO);
    glDeleteVertexArrays(1, &lightCubeVAO);
    glDeleteBuffers(1, &cubeVBO);
    glDeleteBuffers(1, &cubeEBO);

    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();
    return 0;
}

void shaderActivate(Shader& shader)
{
    shader.use();
    shader.setVec3("viewPos", camera.Position);
    shader.setMat4("view", view);
    shader.setMat4("projection", projection);
}


void drawCube(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 model = glm::mat4(1.0f), float r = 1.0f, float g = 1.0f, float b = 1.0f)
{
    lightingShader.use();

    lightingShader.setVec3("material.ambient", glm::vec3(r, g, b));
    lightingShader.setVec3("material.diffuse", glm::vec3(r, g, b));
    lightingShader.setVec3("material.specular", glm::vec3(0.5f, 0.5f, 0.5f));
    lightingShader.setFloat("material.shininess", 32.0f);

    lightingShader.setMat4("model", model);

    glBindVertexArray(cubeVAO);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
}



void floor(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& floor)
{
    

    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(20, .05, 20));
    translate = glm::translate(model, glm::vec3(-0.5, -20.80, -.80));
    model = alTogether * scale * translate;
    lightingShaderWithTexture.setMat4("model", model );
    lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1, 0.753, 0.008));
    glBindVertexArray(cubeVAO);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);



}


void fence(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& floor)
{
    

    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(5, 1.5, 0.1));
    translate = glm::translate(model, glm::vec3(1.0, -0.67, 39));
    model = alTogether * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    floor.drawCubeWithTexture(lightingShaderWithTexture, model);
    //shaderActivate(lightingShader);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(5, 1.5, 0.1));
    translate = glm::translate(model, glm::vec3(-2, -0.67, 39));
    model = alTogether * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    floor.drawCubeWithTexture(lightingShaderWithTexture, model);


    glm::mat4 rotate = glm::mat4(1.0f);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(.1, 1.5, 20.0));
    
    translate = glm::translate(model, glm::vec3(-100, -0.67, -.8));
    model = alTogether *  scale * translate * rotate;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    floor.drawCubeWithTexture(lightingShaderWithTexture, model);



    
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(.1, 1.5, 20.0));

    translate = glm::translate(model, glm::vec3(99.25, -0.67, -.8));
    model = alTogether * scale * translate * rotate;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    floor.drawCubeWithTexture(lightingShaderWithTexture, model);



    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(20, 1.5, 0.1));
    translate = glm::translate(model, glm::vec3(-0.5, -0.67,-160));
    model = alTogether * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    floor.drawCubeWithTexture(lightingShaderWithTexture, model);
 

}






void stall1(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& stall1roof, Cube& stallstick, Cube& stallbase, Cube& stallbase2)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(1.0, 0.5, 0.05));
    //translate = glm::translate(model, glm::vec3(-0.5, 0, 0));
    translate = glm::translate(model, glm::vec3(1.5, 0, 0));
    model = alTogether * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.772, 0.741, 0.486);
    shaderActivate(lightingShaderWithTexture);
    stallbase2.drawCubeWithTexture(lightingShaderWithTexture, model);

    //baseback
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(1.0, 1.74, 0.05));
    //translate = glm::translate(model, glm::vec3(-0.5, 0, -15.0));
    translate = glm::translate(model, glm::vec3(1.5, 0, -15.0));
    model = alTogether * scale * translate;
    stallbase.drawCubeWithTexture(lightingShaderWithTexture, model);

    //baseleft
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(.05, 0.5, 0.5 + .24));
    //translate = glm::translate(model, glm::vec3(-5 - 2 - 2 - 1.0, 0, -1.0));
    translate = glm::translate(model, glm::vec3(-5 - 2 - 2 - 1.0+2+20+5+8+5, 0, -1.0));
    model = alTogether * scale * translate;
    stallbase2.drawCubeWithTexture(lightingShaderWithTexture, model);


    //baseright
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(.05, 0.5, 0.5 + .24));
    //translate = glm::translate(model, glm::vec3(+5 + 2 + 2, 0, -1.0));
    translate = glm::translate(model, glm::vec3(+5 + 2 + 2+2+20+5+8+5, 0, -1.0));
    model = alTogether * scale * translate;
    stallbase2.drawCubeWithTexture(lightingShaderWithTexture, model);

    //stick1
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.1, 1.4, 0.1));
    //translate = glm::translate(model, glm::vec3(-2 - 2 - 1.0, 0, -2 + 1));
    translate = glm::translate(model, glm::vec3(-2 - 2 - 1.0+20, 0, -2 + 1));
    model = alTogether * translate2 * scale * translate;
    stallstick.drawCubeWithTexture(lightingShaderWithTexture, model);


    //stick2
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.1, 1.4, 0.1));
    //translate = glm::translate(model, glm::vec3(2 + 2 + .0, 0, -2 + 1));
    translate = glm::translate(model, glm::vec3(2 + 2 + 2.0+18, 0, -2 + 1));
    model = alTogether * translate2 * scale * translate;
    stallstick.drawCubeWithTexture(lightingShaderWithTexture, model);


    //stick3
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.1, 1.7, 0.1));
    //translate = glm::translate(model, glm::vec3(2 + 2 + .0, 0, -2 - 4 - 1));
    translate = glm::translate(model, glm::vec3(2 + 2 + 2.0+18, 0, -2 - 4 - 1));
    model = alTogether * translate2 * scale * translate;
    stallstick.drawCubeWithTexture(lightingShaderWithTexture, model);


    //stick3
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.1, 1.7, 0.1));
    //translate = glm::translate(model, glm::vec3(-5.0, 0, -2 - 4 - 1));
    translate = glm::translate(model, glm::vec3(-5.0+2+18, 0, -2 - 4 - 1));
    model = alTogether * translate2 * scale * translate;
    stallstick.drawCubeWithTexture(lightingShaderWithTexture, model);

    //roof
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(1.0, 0.05, 1.13));
    //translate = glm::translate(model, glm::vec3(-0.5, 25, -1.2));
    translate = glm::translate(model, glm::vec3(-0.5+2, 25, -1.2));
    float angle = glm::radians(25.0f); // Replace 45.0f with the desired angle in degrees
    glm::mat4 rotateX = glm::rotate(glm::mat4(1.0f), angle, glm::vec3(1.0f, 0.0f, 0.0f));
    model = alTogether * rotateX * scale * translate;
    stall1roof.drawCubeWithTexture(lightingShaderWithTexture, model);

}


void stall2(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& stall1roof, Cube& stallstick, Cube& stallbase, Cube& stallbase2)
{
    //basefront
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(1.0, 0.5, 0.05));
    translate = glm::translate(model, glm::vec3(-0.5, 0, 0));
    model = alTogether * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.501, 0.109, 0.039);
    shaderActivate(lightingShaderWithTexture);
    stall1roof.drawCubeWithTexture(lightingShaderWithTexture, model);

    //baseback
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(1.0, 1.4, 0.05));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -15.0));
    model = alTogether * scale * translate;
    stallbase2.drawCubeWithTexture(lightingShaderWithTexture, model);

    //baseleft
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(.05, 0.5, 0.5 + .24));
    translate = glm::translate(model, glm::vec3(-5 - 2 - 2 - 1.0, 0, -1.0));
    model = alTogether * scale * translate;
    stall1roof.drawCubeWithTexture(lightingShaderWithTexture, model);

    //baseright
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(.05, 0.5, 0.5 + .24));
    translate = glm::translate(model, glm::vec3(+5 + 2 + 2, 0, -1.0));
    model = alTogether * scale * translate;
    stall1roof.drawCubeWithTexture(lightingShaderWithTexture, model);

    //stick1
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.1, 1.4, 0.1));
    translate = glm::translate(model, glm::vec3(-2 - 2 - 1.0, 0, -2 + 1));
    model = alTogether * translate2 * scale * translate;
    stallstick.drawCubeWithTexture(lightingShaderWithTexture, model);


    //stick2
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.1, 1.4, 0.1));
    translate = glm::translate(model, glm::vec3(2 + 2 + .0, 0, -2 + 1));
    model = alTogether * translate2 * scale * translate;
    stallstick.drawCubeWithTexture(lightingShaderWithTexture, model);


    //stick3
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.1, 1.4, 0.1));
    translate = glm::translate(model, glm::vec3(2 + 2 + .0, 0, -2 - 4 - 1));
    model = alTogether * translate2 * scale * translate;
    stallstick.drawCubeWithTexture(lightingShaderWithTexture, model);


    //stick3
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.1, 1.4, 0.1));
    translate = glm::translate(model, glm::vec3(-5.0, 0, -2 - 4 - 1));
    model = alTogether * translate2 * scale * translate;
    stallstick.drawCubeWithTexture(lightingShaderWithTexture, model);

    //roof1
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(1.2, 0.05, 1.1 / 2.0));
    translate = glm::translate(model, glm::vec3(-0.5, 20 + 1, -2.5));
    float angle = glm::radians(40.0f); // Replace 45.0f with the desired angle in degrees
    glm::mat4 rotateX = glm::rotate(glm::mat4(1.0f), angle, glm::vec3(1.0f, 0.0f, 0.0f));
    model = alTogether * rotateX * scale * translate;
    stall1roof.drawCubeWithTexture(lightingShaderWithTexture, model);

    //roof2
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(1.2, 0.05, 1.1 / 2.0));
    translate = glm::translate(model, glm::vec3(-0.5, 20 + 10, .5));
    angle = glm::radians(-40.0f); // Replace 45.0f with the desired angle in degrees
    rotateX = glm::rotate(glm::mat4(1.0f), angle, glm::vec3(1.0f, 0.0f, 0.0f));
    model = alTogether * rotateX * scale * translate;
    stall1roof.drawCubeWithTexture(lightingShaderWithTexture, model);


}



void stall1call(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& stall1roof, Cube& stallstick, Cube& stallbase, Cube& stallbase2)
{
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model, rotateMatrix;
    glm::mat4 identityMatrix = glm::mat4(1.0f);



    translateMatrix = glm::translate(identityMatrix, glm::vec3(-8.0, -1, 7));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(2.5, 2.5, 2.5));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translateMatrix * rotateYMatrix * scaleMatrix;

    lightingShaderWithTexture.setMat4("model", model);
    stall1(cubeVAO, lightingShader, model, lightingShaderWithTexture, stall1roof, stallstick, stallbase, stallbase2);


    translateMatrix = glm::translate(identityMatrix, glm::vec3(-8.0, -1, 4));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(2.5, 2.5, 2.5));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translateMatrix * rotateYMatrix * scaleMatrix;

    lightingShaderWithTexture.setMat4("model", model);
    stall1(cubeVAO, lightingShader, model, lightingShaderWithTexture, stall1roof, stallstick, stallbase, stallbase2);

    translateMatrix = glm::translate(identityMatrix, glm::vec3(-8.0,-1, 1));
    
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(2.5,2.5,2.5));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translateMatrix * rotateYMatrix * scaleMatrix;

    lightingShaderWithTexture.setMat4("model", model);
    stall1(cubeVAO, lightingShader, model, lightingShaderWithTexture, stall1roof, stallstick, stallbase, stallbase2);


    

    translateMatrix = glm::translate(identityMatrix, glm::vec3(-8.0, -1, -2));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(2.5, 2.5, 2.5));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translateMatrix * rotateYMatrix * scaleMatrix;

    lightingShaderWithTexture.setMat4("model", model);
    stall1(cubeVAO, lightingShader, model, lightingShaderWithTexture, stall1roof, stallstick, stallbase, stallbase2);


    translateMatrix = glm::translate(identityMatrix, glm::vec3(-8.0, -1, -5));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(2.5, 2.5, 2.5));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translateMatrix * rotateYMatrix * scaleMatrix;

    lightingShaderWithTexture.setMat4("model", model);
    stall1(cubeVAO, lightingShader, model, lightingShaderWithTexture, stall1roof, stallstick, stallbase, stallbase2);


    translateMatrix = glm::translate(identityMatrix, glm::vec3(-8.0, -1, -8));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(2.5, 2.5, 2.5));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translateMatrix * rotateYMatrix * scaleMatrix;

    lightingShaderWithTexture.setMat4("model", model);
    stall1(cubeVAO, lightingShader, model, lightingShaderWithTexture, stall1roof, stallstick, stallbase, stallbase2);
}



void stall2call(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& stall1roof, Cube& stallstick, Cube& stallbase, Cube& stallbase2)
{
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model, rotateMatrix;
    glm::mat4 identityMatrix = glm::mat4(1.0f);



    translateMatrix = glm::translate(identityMatrix, glm::vec3(8.0, -1, -14));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(2.5, 2.5, 2.5));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(-90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translateMatrix * rotateYMatrix * scaleMatrix;

    lightingShaderWithTexture.setMat4("model", model);
    stall2(cubeVAO, lightingShader, model, lightingShaderWithTexture, stall1roof, stallstick, stallbase, stallbase2);

    translateMatrix = glm::translate(identityMatrix, glm::vec3(8.0, -1, -10));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(2.5, 2.5, 2.5));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(-90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translateMatrix * rotateYMatrix * scaleMatrix;

    lightingShaderWithTexture.setMat4("model", model);
    stall2(cubeVAO, lightingShader, model, lightingShaderWithTexture, stall1roof, stallstick, stallbase, stallbase2);


    translateMatrix = glm::translate(identityMatrix, glm::vec3(8.0, -1, -6));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(2.5, 2.5, 2.5));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(-90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translateMatrix * rotateYMatrix * scaleMatrix;

    lightingShaderWithTexture.setMat4("model", model);
    stall2(cubeVAO, lightingShader, model, lightingShaderWithTexture, stall1roof, stallstick, stallbase, stallbase2);


    translateMatrix = glm::translate(identityMatrix, glm::vec3(8.0, -1, -2));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(2.5, 2.5, 2.5));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(-90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translateMatrix * rotateYMatrix * scaleMatrix;

    lightingShaderWithTexture.setMat4("model", model);
    stall2(cubeVAO, lightingShader, model, lightingShaderWithTexture, stall1roof, stallstick, stallbase, stallbase2);

    translateMatrix = glm::translate(identityMatrix, glm::vec3(8.0, -1, 2));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(2.5, 2.5, 2.5));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(-90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translateMatrix * rotateYMatrix * scaleMatrix;

    lightingShaderWithTexture.setMat4("model", model);
    stall2(cubeVAO, lightingShader, model, lightingShaderWithTexture, stall1roof, stallstick, stallbase, stallbase2);


    
}



void drawCylinder(Shader& lightingShader, glm::vec3 color,glm::mat4 altogether)
{
    lightingShader.use();
    // building model matrix
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, model;
    model = glm::mat4(1.0);
    model = model * altogether;
    lightingShader.setMat4("model", model);
    lightingShader.setVec3("curveColor", color);

    int n = 36;

    vector<float>v;
    vector<int>in;

    for (int i = 0; i < n; i++)
    {
        float ang = (i * 2 * 3.1416) / n;
        float x1 = .5 * cos(ang), z1 = .5 * sin(ang), y1 = .5;
        v.push_back(x1); v.push_back(y1); v.push_back(z1);
        v.push_back(color.r); v.push_back(color.g); v.push_back(color.b);
        float x2 = .5 * cos(ang), z2 = .5 * sin(ang), y2 = -.5;

        v.push_back(x2); v.push_back(y2); v.push_back(z2);
        v.push_back(color.r); v.push_back(color.g); v.push_back(color.b);
    }


    for (int i = 0; i < n - 1; i++)
    {
        in.push_back(i * 2);
        in.push_back((i * 2) + 1);
        in.push_back((i + 1) * 2);
        in.push_back((i + 1) * 2);
        in.push_back((i * 2) + 1);
        in.push_back((i + 1) * 2 + 1);

    }

    in.push_back((n - 1) * 2);
    in.push_back((n - 1) * 2 + 1);
    in.push_back(0);
    in.push_back(0);
    in.push_back((n - 1) * 2 + 1);
    in.push_back(1);



    unsigned int VAO, VBO, EBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, v.size() * sizeof(float), v.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, in.size() * sizeof(float), in.data(), GL_STATIC_DRAW);


    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // vertex normal attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)12);
    glEnableVertexAttribArray(1);

    // setting up materialistic property
    lightingShader.setVec3("material.ambient", color);
    lightingShader.setVec3("material.diffuse", color);
    lightingShader.setVec3("material.specular", glm::vec3(0.5f, 0.5f, 0.5f));
    lightingShader.setFloat("material.shininess", 32.0f);




    glDrawElements(GL_TRIANGLES, in.size(), GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);
    glDeleteVertexArrays(1, &VAO);
}

void canopi(Shader& lightingShader, glm::vec3 color,glm::mat4 altogether) {

    shaderActivate(lightingShader);
    glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model;
   
    model = identityMatrix * altogether;
    lightingShader.setMat4("model", model);

    lightingShader.setVec3("curveColor", glm::vec3(1.0,1.0,1.0));
    can->draw(lightingShader, model);
    lightingShader.setVec3("curveColor", glm::vec3(1.0, 1.0, 1.0));

}

void canopi1(Shader& lightingShader, glm::vec3 color, glm::mat4 altogether) {

    shaderActivate(lightingShader);
    glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model;

    model = identityMatrix * altogether;
    lightingShader.setMat4("model", model);

    lightingShader.setVec3("curveColor", glm::vec3(1.0, 1.0, 1.0));
    can1->draw(lightingShader, model);
    lightingShader.setVec3("curveColor", glm::vec3(1.0, 1.0, 1.0));

}


void cotton(Shader& lightingShader, glm::vec3 color, glm::mat4 altogether) {

    shaderActivate(lightingShader);
    glm::mat4 identityMatrix = glm::mat4(1.0f); 
    glm::mat4 id = glm::mat4(1.0f);

    if (startRotation1)
    {
        rtsp += .5;
        id = glm::rotate(id, glm::radians(rtsp), glm::vec3(0.0f, 1.0f, 0.0f));

    }
    // make sure to initialize matrix to identity matrix first
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model;

    translateMatrix = glm::translate(identityMatrix, glm::vec3(-2.0, 1-.3, 0));
    model = identityMatrix * altogether * translateMatrix *id;
    lightingShader.setMat4("model", model);

    lightingShader.setVec3("curveColor", glm::vec3(1.0, 1.0, 1.0));
    can2->draw(lightingShader, model);
    lightingShader.setVec3("curveColor", glm::vec3(1.0, 1.0, 1.0));


    translateMatrix = glm::translate(identityMatrix, glm::vec3(-2.0, -.2, 0));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1.5, 1.5, 1.5));

    model = altogether * translateMatrix * scaleMatrix;

    // Draw cylinder
    drawCylinder(lightingShader, glm::vec3(0.8902, 0.2941, 0.5010), model);

    lightingShader.setVec3("curveColor", glm::vec3(1.0, 1.0, 1.0));
}

void polls(Shader& lightingShader, glm::vec3 color, glm::mat4 altogether) {

    shaderActivate(lightingShader);
    glm::mat4 identityMatrix = glm::mat4(1.0f);
   
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model;

    

    translateMatrix = glm::translate(identityMatrix, glm::vec3(0,4.8,-5));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(.25, 11.5, .25));

    model = altogether * translateMatrix * scaleMatrix;

    // Draw cylinder
    drawCylinder(lightingShader, glm::vec3(0.8902, 0.2941, 0.5010), model);

    lightingShader.setVec3("curveColor", glm::vec3(1.0, 1.0, 1.0));


    translateMatrix = glm::translate(identityMatrix, glm::vec3(0, 4.8, -5));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(.25, 11.5, .25));

    model = altogether * translateMatrix * scaleMatrix;

    // Draw cylinder
    drawCylinder(lightingShader, glm::vec3(0.8902, 0.2941, 0.5010), model);

    lightingShader.setVec3("curveColor", glm::vec3(1.0, 1.0, 1.0));


}

void canopistall(unsigned int& cubeVAO,Shader& lightingShader, glm::vec3 color, glm::mat4 alTogether, Cube& tex)
{

    shaderActivate(lightingShader);

    glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model;

    translateMatrix = glm::translate(identityMatrix, glm::vec3(0, 3.5, -13.5));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(2.5, 2.0, 1.5));

    model = alTogether * translateMatrix * scaleMatrix;

    lightingShader.setMat4("model", model);

    // Draw canopi
    lightingShader.setVec3("curveColor", color);
    canopi(lightingShader, glm::vec3(1.0f, 0.0f, 0.0f) , model);
    lightingShader.setVec3("curveColor", glm::vec3(1, 1, 1));

    translateMatrix = glm::translate(identityMatrix, glm::vec3(0, 1, -13.5));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(4.5, 4, 4));

    model = alTogether * translateMatrix * scaleMatrix;

    // Draw cylinder
    drawCylinder(lightingShader, glm::vec3(1,1,0),model);

    lightingShader.setVec3("curveColor", glm::vec3(1.0, 1.0, 1.0));


    glm::mat4 model1 = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model1, glm::vec3(1.2, 1.5, 0.1));
    translate = glm::translate(model1, glm::vec3(-0.5, -0.67, -115));
    model1 = alTogether * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShader);
    tex.drawCubeWithTexture(lightingShader, model1);
   


}



void canopistall2(unsigned int& cubeVAO, Shader& lightingShader, glm::vec3 color, glm::mat4 alTogether, Cube& tex)
{

    shaderActivate(lightingShader);

    glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model;

    translateMatrix = glm::translate(identityMatrix, glm::vec3(5, 1.0, -13.5));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1.5, 1.0, 1));

    model = alTogether * translateMatrix * scaleMatrix;

    lightingShader.setMat4("model", model);

    // Draw canopi
    lightingShader.setVec3("curveColor", color);
    canopi(lightingShader, glm::vec3(1.0f, 0.0f, 0.0f), model);
    lightingShader.setVec3("curveColor", glm::vec3(1, 1, 1));

    translateMatrix = glm::translate(identityMatrix, glm::vec3(5, 0, -13.5));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(2.5, 2, 2));

    model = alTogether * translateMatrix * scaleMatrix;

    // Draw cylinder
    drawCylinder(lightingShader, glm::vec3(0, 1, 1), model);

    lightingShader.setVec3("curveColor", glm::vec3(1.0, 1.0, 1.0));


    glm::mat4 model1 = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model1, glm::vec3(.5, 1.0, 0.1));
    translate = glm::translate(model1, glm::vec3(+0.5+5+4, -.98, -125));
    model1 = alTogether * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShader);
    tex.drawCubeWithTexture(lightingShader, model1);

    lightingShader.setVec3("curveColor", glm::vec3(1.0, 1.0, 1.0));

}


void canopistall3(unsigned int& cubeVAO, Shader& lightingShader, glm::vec3 color, glm::mat4 alTogether, Cube& tex)
{

    shaderActivate(lightingShader);

    glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model;

    translateMatrix = glm::translate(identityMatrix, glm::vec3(-5, 1.0, -13.5));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1.5, 1.0, 1));

    model = alTogether * translateMatrix * scaleMatrix;

    lightingShader.setMat4("model", model);

    // Draw canopi
    lightingShader.setVec3("curveColor", color);
    canopi(lightingShader, glm::vec3(1.0f, 0.0f, 0.0f), model);
    lightingShader.setVec3("curveColor", glm::vec3(1, 1, 1));

    translateMatrix = glm::translate(identityMatrix, glm::vec3(-5, 0, -13.5));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(2.5, 2, 2));

    model = alTogether * translateMatrix * scaleMatrix;

    // Draw cylinder
    drawCylinder(lightingShader, glm::vec3(0, 1, 1), model);

    lightingShader.setVec3("curveColor", glm::vec3(1.0, 1.0, 1.0));


    glm::mat4 model1 = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model1, glm::vec3(.5, 1.0, 0.1));
    translate = glm::translate(model1, glm::vec3(-0.5 - 5 - 4 - .5 - .3, -.98, -125));
    model1 = alTogether * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShader);
    tex.drawCubeWithTexture(lightingShader, model1);

    lightingShader.setVec3("curveColor", glm::vec3(1.0, 1.0, 1.0));

}

void drawCircle(GLfloat x, GLfloat y, GLfloat z, GLfloat radius, GLint numberOfSides) {
    // Generate vertex data for the circle
    GLint numberOfVertices = numberOfSides + 1;
    GLfloat doublePi = 2.0f * 3.14159265359f;

    std::vector<GLfloat> vertices;

    for (int i = 0; i < numberOfVertices; i++) {
        vertices.push_back(x + (radius * cos(i * doublePi / numberOfSides))); // X coordinate
        vertices.push_back(z + (radius * sin(i * doublePi / numberOfSides))); // Y coordinate
        vertices.push_back(y);                                               // Z coordinate
    }

    // Create and bind a VAO
    unsigned int VAO, VBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    glBindVertexArray(VAO);

    // Bind and upload vertex data to VBO
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(GLfloat), vertices.data(), GL_STATIC_DRAW);

    // Define vertex attribute layout
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (void*)0);
    glEnableVertexAttribArray(0);

    // Draw the circle as a line loop
    glDrawArrays(GL_LINE_LOOP, 0, numberOfVertices);

    // Cleanup
    glBindVertexArray(0);
    glDeleteBuffers(1, &VBO);
    glDeleteVertexArrays(1, &VAO);
}


void drawCircleflat(GLfloat x, GLfloat y, GLfloat z, GLfloat radius, GLint numberOfSides) {
    // Generate vertex data for the circle
    GLint numberOfVertices = numberOfSides + 1;
    GLfloat doublePi = 2.0f * 3.14159265359f;

    std::vector<GLfloat> vertices;

    for (int i = 0; i < numberOfVertices; i++) {
        
        vertices.push_back(x + (radius * cos(i * doublePi / numberOfSides))); // X coordinate
        vertices.push_back(y); // Y coordinate
        vertices.push_back(z + (radius * sin(i * doublePi / numberOfSides)));                                               // Z coordinate
    }

    // Create and bind a VAO
    unsigned int VAO, VBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    glBindVertexArray(VAO);

    // Bind and upload vertex data to VBO
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(GLfloat), vertices.data(), GL_STATIC_DRAW);

    // Define vertex attribute layout
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (void*)0);
    glEnableVertexAttribArray(0);

    // Draw the circle as a line loop
    glDrawArrays(GL_LINE_LOOP, 0, numberOfVertices);

    // Cleanup
    glBindVertexArray(0);
    glDeleteBuffers(1, &VBO);
    glDeleteVertexArrays(1, &VAO);
}


void drawLine(GLfloat x1, GLfloat y1, GLfloat z1, GLfloat x2, GLfloat y2, GLfloat z2) {
    // Vertex data for the line
    GLfloat vertices[] = {
        x1, y1, z1,
        x2, y2, z2
    };

    // Create and bind a VAO
    unsigned int VAO, VBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    glBindVertexArray(VAO);

    // Bind and upload vertex data to VBO
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    // Define vertex attribute layout
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (void*)0);
    glEnableVertexAttribArray(0);

    // Draw the line
    glDrawArrays(GL_LINES, 0, 2);

    // Cleanup
    glBindVertexArray(0);
    glDeleteBuffers(1, &VBO);
    glDeleteVertexArrays(1, &VAO);
}


void drawFlywheel(glm::vec3 translation, glm::vec3 scale, float rotationAngle, glm::vec3 color, glm::mat4 alTogether, glm::mat4 alTogether2, Shader& lightingShaderWithTexture, unsigned int& cubeVAO) {
    // Apply transformations
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 model1 = glm::mat4(1.0f);
    glm::mat4 model3 = glm::mat4(1.0f);
    glm::mat4 model4 = glm::mat4(1.0f);
    model = glm::translate(model, glm::vec3(-2,-2+1.5+4,-6));
    model = glm::rotate(model, glm::radians(45.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model1 = glm::translate(model1, glm::vec3(2,2+1.5, -2));// Apply translation
    // Apply rotation
    model = glm::scale(model, glm::vec3(3,3,2));
    model1 = glm::scale(model1, glm::vec3(3, 3, 2));

    model = alTogether * model ;
    model1 = alTogether2 * model1 ;

    if (isFlywheelRotating)
    {
        glm::mat4 rotateFlywheelMatrix = glm::rotate(model4, glm::radians(flywheelRotation), glm::vec3(0.0f, 0.0f, 1.0f));
        model3 = model*rotateFlywheelMatrix;
    }
    else
    {
        model3 = model;
    }
    lightingShaderWithTexture.setMat4("model", model3);

    lightingShaderWithTexture.setVec3("curveColor", color);


    
    // Set line width

    glLineWidth(10.0f);

    // Draw circles
    drawCircle(0.0f, -0.4f, 0.0f, 0.5f, 36);  // Small circle behind
    drawCircle(0.0f, 0.0f, 0.0f, 0.5f, 36);   // Large circle in front
    drawCircle(0.0f, -0.4f, 0.0f, 0.9f, 36);  // Larger circle behind
    drawCircle(0.0f, 0.0f, 0.0f, 0.9f, 36);   // Larger circle in front
    drawCircle(0.0f, -0.4f, 0.0f, 0.2f, 36);  // Smallest circle behind
    drawCircle(0.0f, 0.0f, 0.0f, 0.2f, 36);   // Smallest circle in front

    // Draw lines connecting smaller circles
    for (int i = 0; i < 20; i++) {
        GLfloat angleRad = i * 2.0f * 3.1415f / 20;
        GLfloat x1 = 0.2f * cos(angleRad), y1 = 0.2f * sin(angleRad), z1 = -0.4f;
        GLfloat x2 = 0.5f * cos(angleRad), y2 = 0.5f * sin(angleRad), z2 = -0.4f;
        drawLine(x1, y1, z1, x2, y2, z2);
    }
    for (int i = 0; i < 20; i++) {
        GLfloat angleRad = i * 2.0f * 3.1415f / 20;
        GLfloat x1 = 0.2f * cos(angleRad), y1 = 0.2f * sin(angleRad), z1 = 0.0f;
        GLfloat x2 = 0.5f * cos(angleRad), y2 = 0.5f * sin(angleRad), z2 = 0.0f;
        drawLine(x1, y1, z1, x2, y2, z2);
    }

    // Draw lines connecting larger circles
    for (int i = 0; i < 36; i++) {
        GLfloat angleRad = i * 2.0f * 3.1415f / 36;
        GLfloat x1 = 0.5f * cos(angleRad), y1 = 0.5f * sin(angleRad), z1 = -0.4f;
        GLfloat x2 = 0.5f * cos(angleRad), y2 = 0.5f * sin(angleRad), z2 = 0.0f;
        drawLine(x1, y1, z1, x2, y2, z2);
    }

    // Draw lines connecting the largest circles


    // Draw spokes connecting different radii
    for (int i = 0; i < 6; i++) {
        GLfloat angleRad = i * 2.0f * 3.1415f / 6;
        GLfloat x1 = 0.5f * cos(angleRad), y1 = 0.5f * sin(angleRad), z1 = 0.0f;
        GLfloat x2 = 0.9f * cos(angleRad), y2 = 0.9f * sin(angleRad), z2 = 0.0f;
        drawLine(x1, y1, z1, x2, y2, z2);
    }
    for (int i = 0; i < 6; i++) {
        GLfloat angleRad = i * 2.0f * 3.1415f / 6;
        GLfloat x1 = 0.5f * cos(angleRad), y1 = 0.5f * sin(angleRad), z1 = -0.4f;
        GLfloat x2 = 0.9f * cos(angleRad), y2 = 0.9f * sin(angleRad), z2 = -0.4f;

        drawLine(x1, y1, z1, x2, y2, z2);
    }


    for (int i = 0; i < 6; i++) {
        // Calculate the positions of the line's endpoints
        GLfloat angleRad = i * 2.0f * 3.1415f / 6;
        GLfloat x1 = 0.9f * cos(angleRad), y1 = 0.9f * sin(angleRad), z1 = -0.4f;
        GLfloat x2 = 0.9f * cos(angleRad), y2 = 0.9f * sin(angleRad), z2 = 0.0f;

        // Draw the line
        drawLine(x1, y1, z1, x2, y2, z2);

        glm::vec3 midpoint = glm::vec3((x1 + x2) / 2, (y1 + y2) / 2, (z1 + z2) / 2);

    }

    for (int i = 0; i < 6; i++)
    {
        float angle = i * 60.0f; // 60 degrees between each seat
        glm::mat4 seatModel = glm::mat4(1.0f);
        seatModel = glm::translate(seatModel, glm::vec3(0.9f * cos(glm::radians(angle)), 0.9f * sin(glm::radians(angle)), 0.0f)); // move to the correct position on the Ferris wheel
        seatModel = glm::rotate(seatModel, glm::radians(-rotationAngle), glm::vec3(0.0f, 0.0f, 1.0f)); // rotate to match the Ferris wheel's rotation
        seatModel = glm::translate(seatModel, glm::vec3(-0.2f, -0.2f, -0.5f)); // move to the correct height
        seatModel = glm::scale(seatModel, glm::vec3(0.5f, 0.2f, 0.6f));
        lightingShaderWithTexture.setMat4("model", model3 * seatModel);
        lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1, 0, 0));
        glBindVertexArray(cubeVAO);
        glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);
    }

    glm::mat4 seatModel = glm::mat4(1.0f);
    //seatModel = glm::translate(seatModel, glm::vec3(0.9f * cos(glm::radians(angle)), 0.9f * sin(glm::radians(angle)), 0.0f)); // move to the correct position on the Ferris wheel
    //seatModel = glm::rotate(seatModel, glm::radians(-rotationAngle), glm::vec3(0.0f, 0.0f, 1.0f)); // rotate to match the Ferris wheel's rotation
    seatModel = glm::translate(seatModel, glm::vec3(-0.12f, -1.5f, -0.39f)); // move to the correct height
    seatModel = glm::scale(seatModel, glm::vec3(0.25f, 1.5f, 0.15f));
    
    lightingShaderWithTexture.setMat4("model", model * seatModel);
    lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1, 0, 0));
    glBindVertexArray(cubeVAO);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);


   
    //seatModel = glm::translate(seatModel, glm::vec3(0.9f * cos(glm::radians(angle)), 0.9f * sin(glm::radians(angle)), 0.0f)); // move to the correct position on the Ferris wheel
    //seatModel = glm::rotate(seatModel, glm::radians(-rotationAngle), glm::vec3(0.0f, 0.0f, 1.0f)); // rotate to match the Ferris wheel's rotation
    seatModel = glm::translate(seatModel, glm::vec3(-2.2f, 0.0f, -1.5f)); // move to the correct height
    seatModel = glm::scale(seatModel, glm::vec3(5.25f, .05f, 4.5f));

    lightingShaderWithTexture.setMat4("model", model * seatModel);
    lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1, 0, 0));
    glBindVertexArray(cubeVAO);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);


}


void drawFlywheelcall(glm::vec3 translation, glm::vec3 scale, float rotationAngle, glm::vec3 color, glm::mat4 alTogether, glm::mat4 alTogether2, Shader& lightingShaderWithTexture, unsigned int& cubeVAO)
{

    glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model, rotateMatrix;
    if (isFlywheelRotating)
    {

        
        //glm::vec3 flywheelPosition(-4.0f, 2.0f, -4.0f); // Position
        glm::vec3 flywheelPosition(3, 3, 0);
        glm::vec3 flywheelScale(2.0f, 2.0f, 1.0f);     // Scale

        glm::vec3 flywheelColor(0.0f, 0.0f, 10.0f);
        glm::mat4 model1, translateM;
        translateM = glm::mat4(1.0);
        translateM = glm::translate(translateM, glm::vec3(flywheelPosition));
        glm::mat4 rotateFlywheelMatrix = glm::rotate(identityMatrix, glm::radians(flywheelRotation), glm::vec3(0.0f, 0.0f, 1.0f));
        model =alTogether * rotateFlywheelMatrix ;
        model1 =alTogether2 * translateM  ;
        lightingShaderWithTexture.setMat4("model", model);


        drawFlywheel(flywheelPosition, flywheelScale, flywheelRotation, flywheelColor, model, model1, lightingShaderWithTexture, cubeVAO);


        flywheelRotation += 1;
    }

    else
    {
        //glm::vec3 flywheelPosition(-4.0f, 2.0f, -4.0f); // Position
        glm::vec3 flywheelPosition(3, 3, 0);
        glm::vec3 flywheelScale(2.0f, 2.0f, 1.0f);     // Scale

        glm::vec3 flywheelColor(0.0f, 0.0f, 10.0f);

        glm::mat4 model1, translateM;
        translateM = glm::mat4(1.0);
        translateM = glm::translate(translateM, glm::vec3(flywheelPosition));
        model = alTogether2 * translateM;
        model1 =alTogether2 * translateM ;
        lightingShaderWithTexture.setMat4("model", model1);


        drawFlywheel(flywheelPosition, flywheelScale, flywheelRotation, flywheelColor, model, model1, lightingShaderWithTexture, cubeVAO);
    }
}




void drawride(glm::vec3 translation, glm::vec3 scale, float rotationAngle, glm::vec3 color, glm::mat4 alTogether, glm::mat4 alTogether2, Shader& lightingShaderWithTexture, unsigned int& cubeVAO) {
    
    glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model;


    glm::mat4 id1 = glm::mat4(1.0f);
    id1 = glm::translate(id1, glm::vec3(0,-1,0));


    translateMatrix = glm::translate(identityMatrix, glm::vec3(4.9, 2.35, 0));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1.8, .9, 1.5));

    model = alTogether2 * translateMatrix * scaleMatrix ;

    lightingShaderWithTexture.setMat4("model", model);

    // Draw canopi
    lightingShaderWithTexture.setVec3("curveColor", color);
    canopi1(lightingShaderWithTexture, glm::vec3(1.0f, 0.0f, 0.0f), model);
    lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1,1,1));

    glm::mat4 trainModel = glm::mat4(1.0f);
    trainModel = glm::translate(trainModel,glm::vec3(0.0,-3.5,0) ); // move to the correct position on the circular path
    trainModel = glm::scale(trainModel, glm::vec3(0.05f, 5.0f, 0.05f));
    lightingShaderWithTexture.setMat4("model", model * trainModel );
    lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1, 1, 0));
    glBindVertexArray(cubeVAO);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);


    translateMatrix = glm::translate(identityMatrix, glm::vec3(4.9, -0.94, 0));

    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1.5, .2, 1.5));

    model = alTogether2 * translateMatrix * scaleMatrix ;

    lightingShaderWithTexture.setMat4("model", model);
    canopi1(lightingShaderWithTexture, glm::vec3(0.0f, 0.0f, 1.0f), model);
    lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1, 1, 1));
    
    // Draw rail line
    model = glm::mat4(1.0f);
    //model = glm::translate(model, translation); // Apply translation
    // Apply rotation
    model = glm::scale(model, scale);

    model = alTogether2 * model  ;
    lightingShaderWithTexture.setMat4("model", model);

    //lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1,1,0));

    glLineWidth(10.0f);
    drawCircleflat(4.9f, 2.2f, 0.0f, 1.9f, 36);  // Larger circle behind



    static float trainAngle = 0.0f;

    // In your main loop:
    // Increase the train's angle by 0.5 degrees every frame
    if (startRotation) {
        // In your main loop:
        trainAngle += 0.05f; // Increase the train's angle by 0.5 degrees every frame
        if (trainAngle > 360.0f) {
            trainAngle -= 360.0f; // Reset the train's angle when it reaches 360 degrees
        }
    }

    int numSeats = 6; // Number of seats
    float angleIncrement = 360.0f / numSeats; // Angle increment between seats

    for (int i = 0; i < numSeats; i++) {
        glm::mat4 trainModel = glm::mat4(1.0f);
        float seatAngle = trainAngle + (i * angleIncrement); // Calculate the seat angle
        trainModel = glm::translate(trainModel, glm::vec3(4.9+1.7f * cos(glm::radians(seatAngle)), 0.8f, 1.7f * sin(glm::radians(seatAngle)))); // move to the correct position on the circular path
        //trainModel = glm::rotate(trainModel, glm::radians(-seatAngle), glm::vec3(0.0f, 1.0f, 0.0f));
        //trainModel = glm::rotate(trainModel, glm::radians(45.0f), glm::vec3(0.0f, 1.0f, 0.0f));
        trainModel = glm::scale(trainModel, glm::vec3(0.05f, 1.5f, 0.05f));
        lightingShaderWithTexture.setMat4("model", model * trainModel );
        lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1, 1, 0));
        glBindVertexArray(cubeVAO);
        glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);


        //seatAngle = trainAngle + (i * angleIncrement); // Calculate the seat angle
        trainModel = glm::translate(trainModel, glm::vec3(1.7f * cos(glm::radians(seatAngle)), 0.8f, 1.7f * sin(glm::radians(seatAngle)))); // move to the correct position on the circular path
        //trainModel = glm::rotate(trainModel, glm::radians(-seatAngle), glm::vec3(0.0f, 1.0f, 0.0f));
        //trainModel = glm::rotate(trainModel, glm::radians(45.0f), glm::vec3(0.0f, 1.0f, 0.0f));
        trainModel = glm::translate(trainModel, glm::vec3(-8.2f, -1.0f, -3.5f));
        trainModel = glm::scale(trainModel, glm::vec3(18.5f, 0.35f, 12.5f));
        lightingShaderWithTexture.setMat4("model",  model * trainModel );
        lightingShaderWithTexture.setVec3("curveColor", glm::vec3(1, 0, 0));
        glBindVertexArray(cubeVAO);
        glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);

       




    }
}

























void getCurrentTime(int& hours, int& minutes, int& seconds) {
    time_t currentTime = time(nullptr); // Get current UNIX timestamp
    struct tm* timeinfo;
    timeinfo = localtime(&currentTime);
    seconds = timeinfo->tm_sec;
    minutes = timeinfo->tm_min;
    hours = timeinfo->tm_hour;
}
// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
        camera.ProcessKeyboard(FORWARD, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
        camera.ProcessKeyboard(BACKWARD, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
        camera.ProcessKeyboard(LEFT, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
        camera.ProcessKeyboard(RIGHT, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS) {
        camera.ProcessKeyboard(UP, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS) {
        camera.ProcessKeyboard(DOWN, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_T) == GLFW_PRESS) {
        camera.ProcessKeyboard(P_UP, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS) {
        camera.ProcessKeyboard(P_DOWN, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_U) == GLFW_PRESS) {
        camera.ProcessKeyboard(Y_LEFT, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_H) == GLFW_PRESS) {
        camera.ProcessKeyboard(Y_RIGHT, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_I) == GLFW_PRESS) {
        camera.ProcessKeyboard(R_LEFT, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_J) == GLFW_PRESS) {
        camera.ProcessKeyboard(R_RIGHT, deltaTime);
    }


    if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS)
    {
        if (isFlywheelRotating)
        {
            isFlywheelRotating = false;
        }
        else
        {
            isFlywheelRotating = true;
        }
    }
    

    //...
    if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS)
    {
        if (startRotation)
        {
            startRotation = false;
        }
        else
        {
            startRotation = true;
        }
    }

    if (glfwGetKey(window, GLFW_KEY_O) == GLFW_PRESS)
    {
        if (startRotation1)
        {
            startRotation1 = false;
        }
        else
        {
            startRotation1 = true;
        }
    }
 




    /*if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS)
    {
        if (rotateAxis_X) rotateAngle_X -= 0.1;
        else if (rotateAxis_Y) rotateAngle_Y -= 0.1;
        else rotateAngle_Z -= 0.1;
    }*/
    /*if (glfwGetKey(window, GLFW_KEY_I) == GLFW_PRESS) translate_Y += 0.001;*/
    //if (glfwGetKey(window, GLFW_KEY_K) == GLFW_PRESS) translate_Y -= 0.001;
    //if (glfwGetKey(window, GLFW_KEY_L) == GLFW_PRESS) translate_X += 0.001;
    /*if (glfwGetKey(window, GLFW_KEY_J) == GLFW_PRESS) translate_X -= 0.001;*/
    
    //if (glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS) scale_X += 0.001;
    //if (glfwGetKey(window, GLFW_KEY_V) == GLFW_PRESS) scale_X -= 0.001;
    //if (glfwGetKey(window, GLFW_KEY_B) == GLFW_PRESS) scale_Y += 0.001;
    //if (glfwGetKey(window, GLFW_KEY_N) == GLFW_PRESS) scale_Y -= 0.001;
    //if (glfwGetKey(window, GLFW_KEY_M) == GLFW_PRESS) scale_Z += 0.001;
    //if (glfwGetKey(window, GLFW_KEY_U) == GLFW_PRESS) scale_Z -= 0.001;

    if (glfwGetKey(window, GLFW_KEY_X) == GLFW_PRESS)
    {
        rotateAngle_X += 1.0;
        rotateAxis_X = 1.0;
        rotateAxis_Y = 0.0;
        rotateAxis_Z = 0.0;
        angle += 1.0f;
    }
    if (glfwGetKey(window, GLFW_KEY_Y) == GLFW_PRESS)
    {
        rotateAngle_Y += 1.0;
        rotateAxis_X = 0.0;
        rotateAxis_Y = 1.0;
        rotateAxis_Z = 0.0;
        angle += 1.0f;
    }
    if (glfwGetKey(window, GLFW_KEY_Z) == GLFW_PRESS)
    {
        rotateAngle_Z += 1.0;
        rotateAxis_X = 0.0;
        rotateAxis_Y = 0.0;
        rotateAxis_Z = 1.0;
        angle += 1.0f;
    }
    if (glfwGetKey(window, GLFW_KEY_7) == GLFW_PRESS)
    {
        pointlight1.turnAmbientOn();
        pointlight3.turnAmbientOn();
        light1.turnAmbientOn();
        light2.turnAmbientOn();
        light3.turnAmbientOn();

        pointlight1.turnDiffuseOff();
        pointlight3.turnDiffuseOff();
        light1.turnDiffuseOff();
        light2.turnDiffuseOff();
        light3.turnDiffuseOff();

        pointlight1.turnSpecularOff();
        pointlight3.turnSpecularOff();
        light1.turnSpecularOff();
        light2.turnSpecularOff();
        light3.turnSpecularOff();
    }

    if (glfwGetKey(window, GLFW_KEY_8) == GLFW_PRESS)
    {
        pointlight1.turnDiffuseOn();
        pointlight3.turnDiffuseOn();
        light1.turnDiffuseOn();
        light2.turnDiffuseOn();
        light3.turnDiffuseOn();

        pointlight1.turnAmbientOff();
        pointlight3.turnAmbientOff();
        light1.turnAmbientOff();
        light2.turnAmbientOff();
        light3.turnAmbientOff();

        pointlight1.turnSpecularOff();
        pointlight3.turnSpecularOff();
        light1.turnSpecularOff();
        light2.turnSpecularOff();
       light3.turnSpecularOff();

    }

    if (glfwGetKey(window, GLFW_KEY_K) == GLFW_PRESS)
    {
        pointlight1.turnDiffuseOff();
        pointlight3.turnDiffuseOff();
        light1.turnDiffuseOff();
        light2.turnDiffuseOff();
        light3.turnDiffuseOff();

        pointlight1.turnAmbientOff();
        pointlight3.turnAmbientOff();
        light1.turnAmbientOff();
        light2.turnAmbientOff();
        light3.turnAmbientOff();

        pointlight1.turnSpecularOn();
        pointlight3.turnSpecularOn();
        light1.turnSpecularOn();
        light2.turnSpecularOn();
        light3.turnSpecularOn();
    }

    if (glfwGetKey(window, GLFW_KEY_L) == GLFW_PRESS)
    {
        pointlight1.turnDiffuseOn();
        pointlight3.turnDiffuseOn();
        light1.turnDiffuseOn();
        light2.turnDiffuseOn();
        light3.turnDiffuseOn();

        pointlight1.turnAmbientOn();
        pointlight3.turnAmbientOn();
        light1.turnAmbientOn();
        light2.turnAmbientOn();
        light3.turnAmbientOn();

        pointlight1.turnSpecularOn();
        pointlight3.turnSpecularOn();
        light1.turnSpecularOn();
        light2.turnSpecularOn();
        light3.turnSpecularOn();
    }

    if (glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS)
    {
        //pointlight2.turnOn();
        light1.turnOn();
        light2.turnOn();
        light3.turnOn();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_V) == GLFW_PRESS)
    {
        //pointlight2.turnOff();
        light1.turnOff();
        light2.turnOff();
        light3.turnOff();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_B) == GLFW_PRESS)
    {
        pointlight3.turnOn();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_N) == GLFW_PRESS)
    {
        pointlight3.turnOff();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_M) == GLFW_PRESS)
    {
        pointlight1.turnOn();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS)
    {
        pointlight1.turnOff();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }

    if (glfwGetKey(window, GLFW_KEY_5) == GLFW_PRESS)
    {
        /*if (diffuseToggle)
        {*/
        if (pointlight1.isOn())
            pointlight1.turnSpecularOn();
        if (pointlight2.isOn())
            pointlight2.turnSpecularOn();
        if (pointlight3.isOn())
            pointlight3.turnSpecularOn();
        //pointlight4.turnSpecularOn();
        //diffuseToggle = !diffuseToggle;
        //}
    }
    if (glfwGetKey(window, GLFW_KEY_6) == GLFW_PRESS)
    {
        /*if (diffuseToggle)
        {*/
        /*cout << "1 " << pointlight1.isOn() << endl;
        cout << pointlight2.isOn() << endl;
        cout << pointlight3.isOn() << endl;*/
        if (pointlight1.isOn())
            pointlight1.turnSpecularOff();
        if (pointlight2.isOn())
            pointlight2.turnSpecularOff();
        if (pointlight3.isOn())
            pointlight3.turnSpecularOff();
        //pointlight4.turnSpecularOff();
        //diffuseToggle = !diffuseToggle;
        //}
    }
}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}


// glfw: whenever the mouse moves, this callback is called
// -------------------------------------------------------
void mouse_callback(GLFWwindow* window, double xposIn, double yposIn)
{
    float xpos = static_cast<float>(xposIn);
    float ypos = static_cast<float>(yposIn);

    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

    lastX = xpos;
    lastY = ypos;

    camera.ProcessMouseMovement(xoffset, yoffset);
}

// glfw: whenever the mouse scroll wheel scrolls, this callback is called
// ----------------------------------------------------------------------
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.ProcessMouseScroll(static_cast<float>(yoffset));
}

unsigned int loadTexture(char const* path, GLenum textureWrappingModeS, GLenum textureWrappingModeT, GLenum textureFilteringModeMin, GLenum textureFilteringModeMax)
{
    unsigned int textureID;
    glGenTextures(1, &textureID);

    int width, height, nrComponents;
    stbi_set_flip_vertically_on_load(true);
    unsigned char* data = stbi_load(path, &width, &height, &nrComponents, 0);
    if (data)
    {
        GLenum format;
        if (nrComponents == 1)
            format = GL_RED;
        else if (nrComponents == 3)
            format = GL_RGB;
        else if (nrComponents == 4)
            format = GL_RGBA;

        glBindTexture(GL_TEXTURE_2D, textureID);
        glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);
        glGenerateMipmap(GL_TEXTURE_2D);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, textureWrappingModeS);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, textureWrappingModeT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, textureFilteringModeMin);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, textureFilteringModeMax);

        stbi_image_free(data);
    }
    else
    {
        std::cout << "Texture failed to load at path: " << path << std::endl;
        stbi_image_free(data);
    }

    return textureID;
}


