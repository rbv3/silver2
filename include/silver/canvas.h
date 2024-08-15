// Created by Marcos Oliveira <mhco@cin.ufpe.br> on 05/29/2024.
// Copyright (c)
#pragma once

#include <array>
#include <vector>

#include <SFML/Graphics/Color.hpp>
#include <SFML/Graphics/Vertex.hpp>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>

#include "silver/i_node.h"
#include "silver/projection_3d.h"
#include "silver/triangle.h"
#include "silver/window.h"

#include <SFML/Graphics/Texture.hpp>
#include <SFML/Graphics/Image.hpp>

namespace silver {
class Canvas : public INode {
 public:
  glm::vec3 k_{0.2, 0.5, 0.5};
  glm::vec3 lightPos_{10.0, 10.0, 10.0};
  glm::vec3 ambientColor_{255.0, 255.0, 255.0};
  glm::vec3 lightColor_{255.0, 166.0, 246.0};
  
  Canvas(Window* target, Projection3d* projection);

  void AddObject(const std::vector<Triangle>& object, const std::vector<std::vector<std::vector<int>>>& uv);

  void DrawPoint(const glm::vec2& point, const sf::Color& color);
  void DrawPoints(const std::vector<glm::vec2>& points, const sf::Color& color);
  void DrawPoints(const std::vector<glm::vec2>& points,
                  const std::vector<sf::Color>& colors);

  glm::vec3 GetLightIntensity(
    glm::vec3 lightDirection,
    glm::vec3 normal,
    glm::vec3 eyeDirection,
    sf::Color pointColor
  );
  void LoadTextures();

  void Draw();
  void FillTriangle(Triangle& triangle);

  void Clear();

  void OnUpdate(float delta_time) override;

 private:
  Window* target_;
  Projection3d* projection_;
  std::vector<std::vector<Triangle>> objects_;
  std::vector<std::vector<std::vector<std::vector<int>>>> uvs_;

  sf::Texture pernambucoTexture_;
  sf::Texture eyeTexture_;
  sf::Texture lightTexture_;
  sf::Texture normalTexture_;
  sf::Texture paramTexture_;
  sf::Image pernambucoMap_;
  const sf::Uint8* pernambucoMapPtr_;
  sf::Image eyeMap_;
  sf::Image lightMap_;
  sf::Image normalMap_;
  sf::Image paraMap_;
};
}  // namespace silver