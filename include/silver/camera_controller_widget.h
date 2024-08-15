// Created by Marcos Oliveira <mhco@cin.ufpe.br> on 05/25/2024.
// Copyright (c)
#pragma once

#include "silver/camera_3d.h"
#include "silver/canvas.h"
#include "silver/i_widget.h"

namespace silver {

class CameraControllerWidget : public IWidget {
 public:
  explicit CameraControllerWidget(Camera3d* camera, Canvas* canvas);

  void Render() override;

 protected:
  Camera3d* camera_;
  Canvas* canvas_;
};

}  // namespace silver