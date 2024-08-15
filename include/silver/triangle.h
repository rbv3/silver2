// Created by Marcos Oliveira <mhco@cin.ufpe.br> on 06/16/2024.
// Copyright (c)

#pragma once

#include <array>
#include <glm/vec3.hpp>


namespace silver {

struct Triangle {
  std::array<glm::vec3, 3> vertices;
  std::array<glm::vec3, 3> normals;

  Triangle(std::array<glm::vec3, 3> v): vertices(v) {};
};

}  // namespace silver