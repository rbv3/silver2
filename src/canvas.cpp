// Created by Marcos Oliveira <mhco@cin.ufpe.br> on 06/04/2024.
// Copyright (c)

#include <math.h>
#include <algorithm>
#include <array>
#include <format>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>

#include <SFML/Graphics/Color.hpp>
#include <SFML/Graphics/Image.hpp>
#include <SFML/Graphics/Texture.hpp>
#include <SFML/Graphics/Vertex.hpp>
#include <SFML/Graphics/VertexArray.hpp>
#include <glm/geometric.hpp>
#include <glm/gtc/epsilon.hpp>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>

#include "silver/canvas.h"
#include "silver/projection_3d.h"
#include "silver/triangle.h"
#include "silver/window.h"

namespace silver {

namespace internal {
glm::vec3 BaricenterCoordinates(
    const std::array<glm::vec2, 3>& triangule_vertices,
    const glm::vec2& point) {
  glm::vec2 v0 = triangule_vertices[1] - triangule_vertices[0];
  glm::vec2 v1 = triangule_vertices[2] - triangule_vertices[0];
  glm::vec2 v2 = point - triangule_vertices[0];

  float d00 = glm::dot(v0, v0);
  float d01 = glm::dot(v0, v1);
  float d11 = glm::dot(v1, v1);
  float d20 = glm::dot(v2, v0);
  float d21 = glm::dot(v2, v1);
  float denom = d00 * d11 - d01 * d01;
  float w = (d11 * d20 - d01 * d21) / denom;
  float v = (d00 * d21 - d01 * d20) / denom;
  float u = 1.0f - v - w;
  return glm::vec3{u, v, w};
}

glm::vec2 Baricenter(const std::array<glm::vec2, 3>& triangule_vertices) {
  return glm::vec2{(triangule_vertices[0].x + triangule_vertices[1].x +
                    triangule_vertices[2].x) /
                       3.0,
                   (triangule_vertices[0].y + triangule_vertices[1].y +
                    triangule_vertices[2].y) /
                       3.0};
}

glm::vec3 Normal(Triangle triangle) {
  // https://www.khronos.org/opengl/wiki/Calculating_a_Surface_Normal
  glm::vec3 p0 = triangle.vertices[0];
  glm::vec3 p1 = triangle.vertices[1];
  glm::vec3 p2 = triangle.vertices[2];

  glm::vec3 A = p1 - p0;
  glm::vec3 B = p2 - p0;

  return glm::normalize(glm::vec3{
      A.y * B.z - A.z * B.y,
      A.z * B.x - A.x * B.z,
      A.x * B.y - A.y * B.x,
  });
}

glm::vec3 VectorsDirection(glm::vec3 vec1, glm::vec3 vec2) {
  return glm::normalize(glm::vec3{
      vec1.x - vec2.x,
      vec1.y - vec2.y,
      vec1.z - vec2.z,
  });
}

glm::vec3 vectorToColorValue(glm::vec3 v) {
  // v has values from -255 to 255
  // val is [-1, 1] -> +1 -> val is [0, 2]
  // if val is 0 -> 0
  // if val is 1 -> 127.5 (round to 128)
  // if val is 2 -> 255
  float r = std::round(((v.x / 255.0) + 1) * 127.5);
  float g = std::round(((v.y / 255.0) + 1) * 127.5);
  float b = std::round(((v.z / 255.0) + 1) * 127.5);

  return {r, g, b};
}

std::vector<glm::vec2> ScanLine(glm::vec2 p0, glm::vec2 p1, glm::vec2 p2) {
  if (p0.y > p1.y) {
    std::swap(p0, p1);
  }
  if (p1.y > p2.y) {
    std::swap(p1, p2);
  }
  if (p0.y > p1.y) {
    std::swap(p0, p1);
  }

  int x0 = p0.x;
  int y0 = p0.y;
  int x1 = p1.x;
  int y1 = p1.y;
  int x2 = p2.x;
  int y2 = p2.y;

  auto interpolate = [](int y, int y0, int x0, int y1, int x1) {
    return x0 + (x1 - x0) * (y - y0) / (y1 - y0);
  };

  float a01 = (y1 != y0) ? (x1 - x0) / (1.0f * (y1 - y0)) : 0.0f;
  float a02 = (y2 != y0) ? (x2 - x0) / (1.0f * (y2 - y0)) : 0.0f;

  int start_x = x0;
  int end_x = x1;

  std::vector<glm::vec2> points;
  for (int y = y0; y <= y1; y++) {
    int start_x = (y1 != y0) ? interpolate(y, y0, x0, y1, x1) : x1;
    int end_x = (y2 != y0) ? interpolate(y, y0, x0, y2, x2) : x2;
    //start_x += a01;
    //end_x += a02;
    if (start_x > end_x)
      std::swap(start_x, end_x);
    for (int x = start_x; x <= end_x; x++) {
      points.emplace_back(x, y);
    }
  }

  float a12 = (y2 != y1) ? (x2 - x1) / (1.0f * (y2 - y1)) : 0.0f;
  start_x = x1;
  end_x = x2;

  for (int y = y1; y <= y2; y++) {
    int start_x = (y1 != y2) ? interpolate(y, y1, x1, y2, x2) : x2;
    int end_x = (y2 != y0) ? interpolate(y, y0, x0, y2, x2) : x2;
    //start_x += a12;
    //end_x += a02;
    if (start_x > end_x)
      std::swap(start_x, end_x);
    for (int x = start_x; x <= end_x; x++) {
      points.emplace_back(x, y);
    }
  }
  return points;
}
sf::Color uvToTextureColor(int u, int v, int maxU, int maxV,
                           sf::Image& textureMap) {
  int textureWidth = textureMap.getSize().x - 1; // avoid out of range when getPixel [0, size-1]
  int textureHeight = textureMap.getSize().y - 1; // avoid out of range when getPixel [0, size-1]

  int textureX = round((float(u) / maxU) * textureWidth);
  int textureY = round((float(v) / maxV) * textureHeight);

  return textureMap.getPixel(textureX, textureY);
}

int closerTo(glm::vec2 p, glm::vec2 p0, glm::vec2 p1, glm::vec2 p2) {
  float dist0 = (pow(p.x - p0.x, 2) + pow(p.y - p0.y, 2));
  float dist1 = (pow(p.x - p1.x, 2) + pow(p.y - p1.y, 2));
  float dist2 = (pow(p.x - p2.x, 2) + pow(p.y - p2.y, 2));

  bool isBetween0and1 = dist0 <= dist2 && dist1 <= dist2;
  bool isBetween0and2 = dist0 <= dist1 && dist2 <= dist1;
  bool isBetween1and2 = dist1 <= dist0 && dist2 <= dist0;
  if (dist0 >= dist1 && dist0 >= dist2)
    return 0;
  if (dist1 >= dist0 && dist1 >= dist2)
    return 1;
  if (dist2 >= dist1 && dist2 >= dist0)
    return 2;
  return -1;
}
}  // namespace internal

Canvas::Canvas(Window* target, Projection3d* projection)
    : target_(target),
      projection_(projection),
      objects_{},
      uvs_{},
      pernambucoTexture_{},
      eyeTexture_{},
      lightTexture_(),
      normalTexture_(),
      paramTexture_(),
      pernambucoMap_{},
      pernambucoMapPtr_{},
      eyeMap_(),
      lightMap_(),
      normalMap_(),
      paraMap_() {}

void Canvas::LoadTextures() {
  if (!pernambucoTexture_.loadFromFile(
          "/Users/rodrigobeltraovalenca/Documents/silver/assets/textures/"
          "brazil.png")) {
    std::cout << "Failed load" << std::endl;
  } else {
    pernambucoMap_ = pernambucoTexture_.copyToImage();
    pernambucoMapPtr_ = pernambucoMap_.getPixelsPtr();
  }

  int i = 255, j = 255;
}

void Canvas::AddObject(const std::vector<Triangle>& object,
                       const std::vector<std::vector<std::vector<int>>>& uv) {
  LoadTextures();
  objects_.push_back(object);
  uvs_.push_back(uv);
}

void Canvas::DrawPoint(const glm::vec2& point, const sf::Color& color) {
  auto vertex_point =
      sf::Vertex(sf::Vector2f{std::round(point.x), std::round(point.y)}, color);
  target_->window_.draw(&vertex_point, 1, sf::Points);
}

void Canvas::DrawPoints(const std::vector<glm::vec2>& points,
                        const sf::Color& color) {
  sf::VertexArray vertexs(sf::Points, points.size());
  for (int i = 0; i < points.size(); ++i) {
    vertexs[i].position =
        sf::Vector2f{std::round(points[i].x), std::round(points[i].y)};
  }
  target_->window_.draw(vertexs);
}

void Canvas::DrawPoints(const std::vector<glm::vec2>& points,
                        const std::vector<sf::Color>& colors) {
  sf::VertexArray vertexs(sf::Points, points.size());
  for (int i = 0; i < points.size(); ++i) {
    vertexs[i].position =
        sf::Vector2f{std::round(points[i].x), std::round(points[i].y)};
    vertexs[i].color = sf::Color(255, 255, 0);
    vertexs[i].color = sf::Color(colors[i].r, colors[i].g, colors[i].b);
  }
  target_->window_.draw(vertexs);
}

glm::vec3 Canvas::GetLightIntensity(glm::vec3 lightDirection, glm::vec3 normal,
                                    glm::vec3 eyeDirection,
                                    sf::Color pointColor) {
  // I = (Ia * Ka) + (Il * Od * Kd * <N, L> ) + (Il * Ks * <R, V>)
  // componente ambiente = (Ia * Ka)
  // componente difusao = (Il * Od * Kd * <N, L> )
  // componente especular = (Il * Ks * <R, V>)
  // tal que
  // R = 2 * <L, N> - L
  // Ia -> cor da luz ambiente [0, 255]^3
  // Il -> cor da fonte de luz [0, 255]^3
  // Od -> cor do objeto [0,1]^3, atenuador, variados em cada canal, dada pela textura (dada pelo objeto (triangulo baricentro))
  // Ks -> coeficiente especular [0,1]
  // <N, L> vetores q dependem da posiçao do ponto (N vet or normal) (L vetor aponta pra luz)
  // <R, V> vetores q dependem da posiçao do ponto (R vetor reflexao) (V vetor aponta pra observador)
  float ka(k_[0]), kd(k_[1]), ks(k_[2]);
  float sumK = k_[0] + k_[1] + k_[2];

  if (sumK > 1.0) {
    ka /= sumK;
    kd /= sumK;
    ks /= sumK;
  }
  glm::vec3 normalizedNormal = glm::normalize(normal);
  glm::vec3 normalizedLight = glm::normalize(lightDirection);
  glm::vec3 normalizedEye = glm::normalize(eyeDirection);

  float dotLightNormal = glm::dot(normalizedLight, normalizedNormal);

  glm::vec3 lightNormalCross =
      glm::normalize(glm::cross(lightDirection, normal));
  glm::vec3 R = glm::normalize(glm::vec3{
      (2.0 * dotLightNormal) * normalizedNormal.x - normalizedLight.x,
      (2.0 * dotLightNormal) * normalizedNormal.y - normalizedLight.y,
      (2.0 * dotLightNormal) * normalizedNormal.z - normalizedLight.z,
  });

  float dotEyeR = glm::dot(normalizedEye, R);

  sf::Color lightColor = sf::Color(lightColor_.x, lightColor_.y, lightColor_.z);
  sf::Color ambientColor =
      sf::Color(ambientColor_.x, ambientColor_.y, ambientColor_.z);
  glm::vec3 normalizedAmbientLight = {
      (int(ambientColor.r) / 255.0),
      (int(ambientColor.g) / 255.0),
      (int(ambientColor.b) / 255.0),
  };
  glm::vec3 normalizedDirectionalLight = {
      (int(lightColor.r) / 255.0),
      (int(lightColor.g) / 255.0),
      (int(lightColor.b) / 255.0),
  };
  glm::vec3 ambientIntensity = {
      //[0, 255]
      normalizedAmbientLight.x * ka,
      normalizedAmbientLight.y * ka,
      normalizedAmbientLight.z * ka,
  };
  glm::vec3 difusionIntensity = {
      //[-255, 255]
      normalizedDirectionalLight.x * (int(pointColor.r) / 255.0) * kd *
          dotLightNormal,
      normalizedDirectionalLight.y * (int(pointColor.g) / 255.0) * kd *
          dotLightNormal,
      normalizedDirectionalLight.z * (int(pointColor.b) / 255.0) * kd *
          dotLightNormal,
  };
  glm::vec3 specularIntensity = {
      //[-255, 255]
      normalizedDirectionalLight.x * ks * dotEyeR,
      normalizedDirectionalLight.y * ks * dotEyeR,
      normalizedDirectionalLight.z * ks * dotEyeR,
  };

  if (dotLightNormal <= 0)
    difusionIntensity = {0, 0, 0};
  if (dotEyeR <= 0)
    specularIntensity = {0, 0, 0};

  glm::vec3 colorIntensity =
      ambientIntensity + difusionIntensity + specularIntensity;


  // return {
  //     normalizedNormal.x * 255.0,
  //     normalizedNormal.y * 255.0,
  //     normalizedNormal.z * 255.0,
  // };
  return {
      colorIntensity.x * 255.0,
      colorIntensity.y * 255.0,
      colorIntensity.z * 255.0,
  };
}

void Canvas::FillTriangle(Triangle& triangle) {
  std::array<glm::vec2, 3> projections = {
      projection_->Project(triangle.vertices[0]).value(),
      projection_->Project(triangle.vertices[1]).value(),
      projection_->Project(triangle.vertices[2]).value()};

  std::vector<glm::vec2> points =
      internal::ScanLine(projections[0], projections[1], projections[2]);

  constexpr int kNumWhiteTones = 5;
  constexpr int kWhiteTonesDelta = 255 / kNumWhiteTones;

  std::vector<sf::Color> white_tones(kNumWhiteTones);
  for (int i = 0, j = 0; i <= 255 && j < kNumWhiteTones;
       i += kWhiteTonesDelta, ++j) {
    white_tones[j] = sf::Color(i, i, i, 255);
  }

  auto baricenter = internal::Baricenter(projections);

  auto tone_index =
      static_cast<int>(((target_->window_.getSize().y - baricenter.y) /
                        target_->window_.getSize().y) *
                       kNumWhiteTones);

  DrawPoints(points, white_tones[tone_index]);
}

void Canvas::Draw() {
  auto to_index_buffer = [](const glm::vec2& position, int columns) -> size_t {
    size_t x_idx = static_cast<size_t>(std::round(position.x));
    size_t y_idx = static_cast<size_t>(std::round(position.y));
    return x_idx * columns + y_idx;
  };

  std::unordered_map<size_t, float> z_buffer;
  std::unordered_map<size_t, sf::Color> colors;

  auto color_by_point = [](const float& normalized_height) -> sf::Color {
    glm::uint8 gray_scale = static_cast<glm::uint8>(normalized_height * 255.0f);
    return sf::Color{gray_scale, gray_scale, gray_scale, 255};
  };

  std::vector<glm::vec2> render_points;
  std::vector<sf::Color> render_colors;
  auto size = target_->window_.getSize();

  std::map<std::vector<float>, std::vector<int>> pointsToTriangle;
  std::unordered_map<int, glm::vec3> triangleNormals;

  for (auto& object : objects_) {
    int triangleIndex = 0;
    for (auto& triangule : object) {
      auto p0 = projection_->Project(triangule.vertices[0]);
      auto p1 = projection_->Project(triangule.vertices[1]);
      auto p2 = projection_->Project(triangule.vertices[2]);
      auto normal = internal::Normal(triangule);
      triangleNormals.insert({triangleIndex, normal});

      if (p0.has_value() && p1.has_value() && p2.has_value()) {
        std::vector<glm::vec2> points =
            internal::ScanLine(p0.value(), p1.value(), p2.value());

        for (const auto& point : points) {
          glm::vec3 coeff = internal::BaricenterCoordinates(
              {p0.value(), p1.value(), p2.value()}, point);

          glm::vec3 view_point =
              projection_->camera_->Translate(triangule.vertices[0]) * coeff.x +
              projection_->camera_->Translate(triangule.vertices[1]) * coeff.y +
              projection_->camera_->Translate(triangule.vertices[2]) * coeff.z;
          std::vector<float> vec_view_point = {view_point.x, view_point.y,
                                               view_point.z};
          pointsToTriangle[vec_view_point].push_back(triangleIndex);
        }

        triangleIndex++;
      }
    }
  }

  int uvIndex = 0;
  for (auto& object : objects_) {
    int triangleIndex = 0;

    int uvSize = uvs_[uvIndex].size();
    auto lastUvObject = uvs_[uvIndex][uvSize - 1][1];  //i+1, j+1 pair
    auto maxU = lastUvObject[0];
    auto maxV = lastUvObject[1];
    
    for (auto& triangule : object) {
      //each triangle vertice has a unique uv coord
      auto u0 = uvs_[uvIndex][triangleIndex][0][0];
      auto v0 = uvs_[uvIndex][triangleIndex][0][1];
      auto u1 = uvs_[uvIndex][triangleIndex][1][0];
      auto v1 = uvs_[uvIndex][triangleIndex][1][1];
      auto u2 = uvs_[uvIndex][triangleIndex][2][0];
      auto v2 = uvs_[uvIndex][triangleIndex][2][1];

      auto color0 =
          internal::uvToTextureColor(u0, v0, maxU, maxV, pernambucoMap_);
      auto color1 =
          internal::uvToTextureColor(u1, v1, maxU, maxV, pernambucoMap_);
      auto color2 =
          internal::uvToTextureColor(u2, v2, maxU, maxV, pernambucoMap_);
      
      // std::cout<<int(color0.r)<<" "<<int(color0.g)<<" "<<int(color0.b)<<" "<<std::endl;
      auto p0 = projection_->Project(triangule.vertices[0]);
      auto p1 = projection_->Project(triangule.vertices[1]);
      auto p2 = projection_->Project(triangule.vertices[2]);

      auto normal = internal::Normal(triangule);

      if (p0.has_value() && p1.has_value() && p2.has_value()) {
        std::vector<glm::vec2> points =
            internal::ScanLine(p0.value(), p1.value(), p2.value());
        int pointIndex = 0;
        for (const auto& point : points) {
          int closerTo = internal::closerTo(point, p0.value(), p1.value(),
                                                 p2.value());
          pointIndex++;
          sf::Color pointTextureColor =
              closerTo == 0 ? color0 : (closerTo == 1 ? color1 : color2);
          // sf::Color pointTextureColor = color0;


          glm::vec3 coeff = internal::BaricenterCoordinates(
              {p0.value(), p1.value(), p2.value()}, point);
          glm::vec3 realWorldCoords = triangule.vertices[0] * coeff.x +
                                      triangule.vertices[1] * coeff.y +
                                      triangule.vertices[2] * coeff.z;

          glm::vec3 view_point =
              projection_->camera_->Translate(triangule.vertices[0]) * coeff.x +
              projection_->camera_->Translate(triangule.vertices[1]) * coeff.y +
              projection_->camera_->Translate(triangule.vertices[2]) * coeff.z;
          size_t matriz_index = to_index_buffer(point, size.y);

          std::vector<float> vec_view_point = {view_point.x, view_point.y,
                                               view_point.z};

          // get avg normal
          auto currentNormal = normal;
          if (!z_buffer.contains(matriz_index) ||
              view_point.z < z_buffer[matriz_index]) {
            // if (pointsToTriangle[vec_view_point].size() > 1) {
            //   glm::vec3 tempNormal = {0, 0, 0};
            //   for (auto index : pointsToTriangle[vec_view_point]) {
            //     auto indexNormal = triangleNormals[index];
            //     tempNormal += indexNormal;
            //   }
            //   tempNormal /= pointsToTriangle[vec_view_point].size();
            //   currentNormal = tempNormal;
            // }

            auto baricenter =
                internal::Baricenter({p0.value(), p1.value(), p2.value()});

            // calculate directions for light
            auto lightDirection =
                internal::VectorsDirection(lightPos_, realWorldCoords);
            auto eyeDirection = internal::VectorsDirection(
                realWorldCoords, projection_->camera_->GetWorldCoordinate());

            auto pointColor = color_by_point((size.y - baricenter.y) / size.y);

            auto colorIntensity = GetLightIntensity(
                lightDirection, currentNormal, eyeDirection, pointTextureColor);
            // std::cout<<colorIntensity.x<<" "<<colorIntensity.y<<" "<<colorIntensity.z<<std::endl;
            colors[matriz_index] = {
                (uint8_t)round(colorIntensity.x),
                (uint8_t)round(colorIntensity.y),
                (uint8_t)round(colorIntensity.z),
            };
            // colors[matriz_index] = pointTextureColor;

            z_buffer[matriz_index] = view_point.z;
          }

          render_points.push_back(point);
        }
      }
      triangleIndex++;
    }
    uvIndex++;
  }

  for (const auto& point : render_points) {
    size_t matriz_index = to_index_buffer(point, size.y);
    render_colors.push_back(colors[matriz_index]);
  }

  DrawPoints(render_points, render_colors);
}

void Canvas::Clear() {
  objects_.clear();
}

void Canvas::OnUpdate(float delta_time) {
  Draw();
}

}  // namespace silver