// Created by Marcos Oliveira <mhco@cin.ufpe.br> on 05/25/2024.
// Copyright (c)

#include <format>

#include <imgui.h>
#include <glm/gtc/type_ptr.hpp>

#include "silver/camera_controller_widget.h"

namespace silver {
CameraControllerWidget::CameraControllerWidget(Camera3d* camera, Canvas* canvas)
    : camera_(camera), canvas_(canvas) {}

void CameraControllerWidget::Render() {
  ImGui::Begin("Camera Controller");
  ImGui::Text("x");
  ImGui::SameLine();
  ImGui::Text("y");
  ImGui::SameLine();
  ImGui::Text("z");
  bool value_changed = false;
  ImGui::BeginDisabled();
  ImGui::SliderFloat3("Cross Normal x Up vector", glm::value_ptr(camera_->u_),
                      0.0, 1.0, "%.2f", ImGuiSliderFlags_NoInput);
  ImGui::EndDisabled();
  value_changed = value_changed ||
                  ImGui::SliderFloat3("Up vector", glm::value_ptr(camera_->v_),
                                      -1.0f, 1.0f, "%.2f");
  value_changed = value_changed ||
                  ImGui::SliderFloat3("Normal", glm::value_ptr(camera_->n_),
                                      -1.0f, 1.0f, "%.2f");
  value_changed = value_changed ||
                  ImGui::SliderFloat3("Focus", glm::value_ptr(camera_->focus_),
                                      -10.0f, 10.0f, "%.2f");
  value_changed =
      value_changed || ImGui::SliderFloat("Distance from projection",
                                          &camera_->distance_from_projection_,
                                          0.1f, 100.00f, "%.2f");
  value_changed =
      value_changed ||
      ImGui::SliderFloat2("hx | hy", glm::value_ptr(camera_->visible_area_),
                          0.5f, 10.f, "%.2f");


  ImGui::SliderFloat3("Light Pos", glm::value_ptr(canvas_->lightPos_),
                      -10.0, 10.0, "%.2f", ImGuiSliderFlags_NoInput);
  ImGui::SliderFloat3("Ka, Kd, Ks", glm::value_ptr(canvas_->k_),
                      0.0, 1.0, "%.2f", ImGuiSliderFlags_NoInput);
  ImGui::SliderFloat3("AmbientColor", glm::value_ptr(canvas_->ambientColor_),
                      0.0, 255.0, "%.2f", ImGuiSliderFlags_NoInput);
  ImGui::SliderFloat3("LightColor", glm::value_ptr(canvas_->lightColor_),
                      0.0, 255.0, "%.2f", ImGuiSliderFlags_NoInput);
                      
  if (value_changed) {
    camera_->Orthonomalize();
  }
  if (ImGui::Button("Reset")) {
    camera_->Reset();
  }
  ImGui::End();
}

}  // namespace silver