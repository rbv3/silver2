// Created by Marcos Oliveira <mhco@cin.ufpe.br> on 05/22/2024.
// Copyright (c)
#pragma once

#include <glm/mat3x3.hpp>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>

namespace silver {

class CameraControllerWidget;
class Projection3d;

struct CameraConfig {
  glm::vec3 focus;
  glm::vec3 n;
  glm::vec3 v;
  glm::vec2 visible_area;
  float distance_from_projection;
};

class Camera3d {
 public:
  explicit Camera3d();

  void Orthonomalize();

  void Reset();
  
  void Set(const CameraConfig& config);

  const glm::vec3 & GetWorldCoordinate();

  glm::vec3 Translate(const glm::vec3& point);

 protected:
  glm::vec3 focus_;
  glm::vec3 n_;
  glm::vec3 v_;
  glm::vec3 u_;
  glm::vec2 visible_area_;
  float distance_from_projection_;

  friend class CameraControllerWidget;
  friend class Projection3d;
};

}  // namespace silver