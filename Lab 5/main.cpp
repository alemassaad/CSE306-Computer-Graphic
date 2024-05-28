#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <cstdio>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

class Vector {
public:
    Vector(double x = 0, double y = 0, double z = 0) : data{x, y, z} {}

    double length_squared() const {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }

    double length() const {
        return std::sqrt(length_squared());
    }

    void normalize() {
        double len = length();
        data[0] /= len;
        data[1] /= len;
        data[2] /= len;
    }

    double operator[](int idx) const { return data[idx]; }
    double& operator[](int idx) { return data[idx]; }

private:
    double data[3];
};

Vector operator+(const Vector &a, const Vector &b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

Vector operator-(const Vector &a, const Vector &b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

Vector operator*(double scalar, const Vector &v) {
    return Vector(scalar * v[0], scalar * v[1], scalar * v[2]);
}

Vector operator*(const Vector &v, double scalar) {
    return Vector(v[0] * scalar, v[1] * scalar, v[2] * scalar);
}

double dot_product(const Vector &a, const Vector &b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vector generate_random_vector() {
    double r1 = static_cast<double>(rand()) / RAND_MAX;
    double r2 = static_cast<double>(rand()) / RAND_MAX;

    double x = cos(2 * M_PI * r1) * std::sqrt(r2 * (1 - r2));
    double y = sin(2 * M_PI * r1) * std::sqrt(r2 * (1 - r2));
    double z = 1 - 2 * r2;
    return Vector(x, y, z);
}

unsigned char* resize_image(const unsigned char* img, int old_w, int old_h, int channels, int new_w, int new_h) {
    unsigned char* resized_img = new unsigned char[new_w * new_h * channels];
    for (int y = 0; y < new_h; ++y) {
        for (int x = 0; x < new_w; ++x) {
            float gx = x * (old_w - 1) / static_cast<float>(new_w - 1);
            float gy = y * (old_h - 1) / static_cast<float>(new_h - 1);
            int gxi = static_cast<int>(gx);
            int gyi = static_cast<int>(gy);

            for (int c = 0; c < channels; ++c) {
                float c00 = img[(gyi * old_w + gxi) * channels + c];
                float c10 = img[(gyi * old_w + gxi + 1) * channels + c];
                float c01 = img[((gyi + 1) * old_w + gxi) * channels + c];
                float c11 = img[((gyi + 1) * old_w + gxi + 1) * channels + c];
                resized_img[(y * new_w + x) * channels + c] = static_cast<unsigned char>(
                    c00 * (1 - (gx - gxi)) * (1 - (gy - gyi)) +
                    c10 * (gx - gxi) * (1 - (gy - gyi)) +
                    c01 * (1 - (gx - gxi)) * (gy - gyi) +
                    c11 * (gx - gxi) * (gy - gyi)
                );
            }
        }
    }
    return resized_img;
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        printf("Invalid arguments.\n");
        return 1;
    }

    const char *input_path = argv[1];
    const char *model_path = argv[2];
    int iterations = atoi(argv[3]);

    int input_w, input_h, input_channels;
    unsigned char *input_image = stbi_load(input_path, &input_w, &input_h, &input_channels, 0);
    if (!input_image) {
        printf("Failed to load input image: %s\n", input_path);
        return 1;
    }

    int model_w, model_h, model_channels;
    unsigned char *model_image = stbi_load(model_path, &model_w, &model_h, &model_channels, 0);
    if (!model_image) {
        printf("Failed to load model image: %s\n", model_path);
        stbi_image_free(input_image);
        return 1;
    }

    if (input_channels != model_channels) {
        printf("Image channels do not match.\n");
        stbi_image_free(input_image);
        stbi_image_free(model_image);
        return 1;
    }

    int new_w = std::min(input_w, model_w);
    int new_h = std::min(input_h, model_h);

    unsigned char *resized_input = resize_image(input_image, input_w, input_h, input_channels, new_w, new_h);
    unsigned char *resized_model = resize_image(model_image, model_w, model_h, model_channels, new_w, new_h);

    size_t total_pixels = new_w * new_h;
    std::vector<std::pair<double, int>> projected_input(total_pixels);
    std::vector<std::pair<double, int>> projected_model(total_pixels);

    Vector pixel, model_pixel, rand_vec;

    for (int iter = 0; iter < iterations; ++iter) {
        rand_vec = generate_random_vector();

        for (size_t i = 0; i < total_pixels; ++i) {
            unsigned char *in_pixel = resized_input + input_channels * i;
            unsigned char *mod_pixel = resized_model + model_channels * i;

            pixel = Vector(in_pixel[0], in_pixel[1], in_pixel[2]);
            model_pixel = Vector(mod_pixel[0], mod_pixel[1], mod_pixel[2]);

            projected_input[i] = {dot_product(pixel, rand_vec), static_cast<int>(i)};
            projected_model[i] = {dot_product(model_pixel, rand_vec), static_cast<int>(i)};
        }

        std::sort(projected_input.begin(), projected_input.end());
        std::sort(projected_model.begin(), projected_model.end());

        for (size_t i = 0; i < total_pixels; ++i) {
            int idx = projected_input[i].second;
            unsigned char *in_pixel = resized_input + input_channels * idx;
            pixel = Vector(in_pixel[0], in_pixel[1], in_pixel[2]);
            Vector adjustment = (projected_model[i].first - projected_input[i].first) * rand_vec;
            pixel = pixel + adjustment;
            in_pixel[0] = static_cast<unsigned char>(pixel[0]);
            in_pixel[1] = static_cast<unsigned char>(pixel[1]);
            in_pixel[2] = static_cast<unsigned char>(pixel[2]);
        }
    }

    stbi_write_png("output.png", new_w, new_h, input_channels, resized_input, 0);

    stbi_image_free(input_image);
    stbi_image_free(model_image);
    delete[] resized_input;
    delete[] resized_model;

    return 0;
}
