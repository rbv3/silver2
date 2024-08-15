// Created by Marcos Oliveira <mhco@cin.ufpe.br> on 05/25/2024.
// Copyright (c)

#include <format>

#include <imgui.h>
#include <glm/gtc/type_ptr.hpp>

#include "silver/illumination_controller_widget.h"

namespace silver {
IlluminationControllerWidget::IlluminationControllerWidget(Canvas* canvas)
    : canvas_(canvas) {}

void IlluminationControllerWidget::Render() {
  ImGui::Begin("Illumination Controller");
  
  ImGui::SliderFloat("Ka", glm::value_ptr(canvas_->ka_),
                      0.0, 1.0, "%.2f", ImGuiSliderFlags_NoInput);
  ImGui::End();
}

}  // namespace silver