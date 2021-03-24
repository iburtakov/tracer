#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "geometry.h"

struct Light {
    Light(const Vec3f &p, const float &i) : position(p), intensity(i) {}
    Vec3f position;
    float intensity;
};


struct Material {
    Material(const Vec2f &a, const Vec3f &color, const float &spec) : albedo(a), diffuse_color(color), specular_exponent(spec) {}
    Material() : albedo(1,0), diffuse_color(), specular_exponent() {}
    Vec2f albedo;
    Vec3f diffuse_color;
    float specular_exponent;
};

// struct Sphere {
//     Vec3f center;
//     float radius;
//     Material material;

//     Sphere(const Vec3f &c, const float &r, const Material &m) : center(c), radius(r), material(m) {}

//     bool ray_intersect(const Vec3f &orig, const Vec3f &dir, float &t0) const {
//         Vec3f L = center - orig;
//         float tca = L*dir;
//         float d2 = L*L - tca*tca;
//         if (d2 > radius*radius) return false;
//         float thc = sqrtf(radius*radius - d2);
//         t0       = tca - thc;
//         float t1 = tca + thc;
//         if (t0 < 0) t0 = t1;
//         if (t0 < 0) return false;
//         return true;
//     }
// };

Vec3f reflect(const Vec3f &I, const Vec3f &N) {
    return I - N*2.f*(I*N);
}


struct Triangle{
    const float EPSILON = 1e-8;
    Vec3f vert0;
    Vec3f vert1;
    Vec3f vert2;
    Material material;

    Triangle(const Vec3f &v0, const Vec3f &v1, const Vec3f &v2, const Material &m) : vert0(v0), vert1(v1), vert2(v2), material(m) {}

    bool ray_intersect(const Vec3f &orig, const Vec3f &dir, float &t) const {
        float u; float v; float w;
        Vec3f edge1 = vert1 - vert0;
        Vec3f edge2 = vert2 - vert0;

        Vec3f pvec = cross(dir, edge2);

        float det = edge1 * pvec;

        if (det < EPSILON)
            return false;
        float inv_det = 1.0f / det;

        Vec3f tvec = orig - vert0;

        u = (tvec * pvec) * inv_det;

        if (u < 0.0 || u > 1.0f)
            return false;

        Vec3f qvec = cross(tvec, edge1);

        v = (dir * qvec) * inv_det;
        if (v < 0.0 || u + v > 1.0f)
            return false;

        t = (edge2 * qvec) * inv_det;

        return true;
    }
};

bool scene_intersect(const Vec3f &orig, const Vec3f &dir, const std::vector<Triangle> &triangles, Vec3f &hit, Vec3f &N, Material &material) {
    float triangles_dist = std::numeric_limits<float>::max();
    for (size_t i=0; i < triangles.size(); i++) {
        float dist_i;
        if (triangles[i].ray_intersect(orig, dir, dist_i) && dist_i < triangles_dist) {
            triangles_dist = dist_i;
            hit = orig + dir*dist_i;
            //N = (hit - spheres[i].center).normalize();
            N = (cross(triangles[i].vert1 - triangles[i].vert0, triangles[i].vert2 - triangles[i].vert0)).normalize();
            // if (N * hit < 0)
            //     N = N * (-1);
            material = triangles[i].material;
        }
    }
    return triangles_dist<1000;
}

Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, const std::vector<Triangle> &triangles, const std::vector<Light> &lights){
    Vec3f point, N;
    Material material;

    if (!scene_intersect(orig, dir, triangles, point, N, material)){
        return Vec3f(0.2, 0.7, 0.8); // background color
    }
    //return material.diffuse_color;
    float diffuse_light_intensity = 0, specular_light_intensity = 0;
    for (size_t i=0; i<lights.size(); i++) {
        Vec3f light_dir      = (lights[i].position - point).normalize();
        float light_distance = (lights[i].position - point).norm();
        //shadows
        Vec3f shadow_orig = light_dir*N < 0 ? point - N*1e-3 : point + N*1e-3; // checking if the point lies in the shadow of the lights[i]
        Vec3f shadow_pt, shadow_N;
        Material tmpmaterial;
        if (scene_intersect(shadow_orig, light_dir, triangles, shadow_pt, shadow_N, tmpmaterial) && (shadow_pt-shadow_orig).norm() < light_distance)
            continue;
        
        diffuse_light_intensity  += lights[i].intensity * std::max(0.f, light_dir*N);
        specular_light_intensity += powf(std::max(0.f, -reflect(-light_dir, N)*dir), material.specular_exponent)*lights[i].intensity;
    }
    return material.diffuse_color * diffuse_light_intensity * material.albedo[0] + Vec3f(1., 1., 1.)*specular_light_intensity * material.albedo[1];
}
/* old cast_ray func
Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, const Triangle &triangle) {
    float triangle_dist = std::numeric_limits<float>::max();
    float t = triangle_dist; //t = dist_i
    if (!triangle.ray_intersect(orig, dir, t)) {
        return Vec3f(0.2, 0.7, 0.8); // background color
    }

    return Vec3f(0.4, 0.4, 0.3);
}
*/
void render(const std::vector<Triangle> &triangles, const std::vector<Light> &lights) {
    const int width    = 1024;
    const int height   = 768;
    const int fov      = M_PI/2.;
    std::vector<Vec3f> framebuffer(width*height);

    #pragma omp parallel for
    for (size_t j = 0; j<height; j++) {
        for (size_t i = 0; i<width; i++) {
            float x =  (2*(i + 0.5)/(float)width  - 1)*tan(fov/2.)*width/(float)height;
            float y = -(2*(j + 0.5)/(float)height - 1)*tan(fov/2.);
            Vec3f dir = Vec3f(x, y, -1).normalize();
            framebuffer[i+j*width] = cast_ray(Vec3f(0,0,0), dir, triangles, lights);
        }
    }

    std::ofstream ofs; // save the framebuffer to file
    ofs.open("./out.ppm");
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height*width; ++i) {
        Vec3f &c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max>1) c = c*(1./max);
        for (size_t j = 0; j<3; j++) {
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    ofs.close();
}

int main() {
    Material      ivory(Vec2f(0.6,  0.3), Vec3f(0.4, 0.4, 0.3),   50.);
    Material red_rubber(Vec2f(0.9,  0.1), Vec3f(0.3, 0.1, 0.1),   50.);
    //Sphere sphere(Vec3f(-3, 0, -16), 2);

    std::vector<Triangle> triangles;
    
    triangles.push_back(Triangle(Vec3f(-3, -3, -15), Vec3f(3, -3, -10), Vec3f(0, 3, -7), red_rubber));
    triangles.push_back(Triangle(Vec3f(-1, -3, -9.6), Vec3f(1, -3, -9.6), Vec3f(0, -1, -9.6), ivory));
    std::vector<Light>  lights;
    lights.push_back(Light(Vec3f(0, 0, -3), 4));

    render(triangles, lights);

    return 0;
}