// Created by Marcos Oliveira <mhco@cin.ufpe.br> on 05/25/2024.
// Copyright (c)
#pragma once

#include "silver/canvas.h"
#include "silver/i_widget.h"

namespace silver {

class IlluminationControllerWidget : public IWidget {
 public:
  explicit IlluminationControllerWidget(Canvas* canvas_);

  void Render() override;

 protected:
  Canvas* canvas_;
};

}  // namespace silver