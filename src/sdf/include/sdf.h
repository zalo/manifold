// Copyright 2023 The Manifold Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#include <functional>

#include "public.h"

namespace manifold {
Mesh LevelSet(std::function<float(glm::vec3)> sdf, Box bounds, float edgeLength,
              float level = 0, bool canParallel = true);
Mesh LevelSetBatch(
    std::function<std::vector<float>(std::vector<glm::vec3>)> sdf, Box bounds,
    float edgeLength, float level = 0, bool canParallel = true);
}  // namespace manifold
