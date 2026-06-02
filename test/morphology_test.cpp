// Copyright 2024 The Manifold Authors.
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

#include "manifold/manifold.h"
#include "test.h"

using namespace manifold;

namespace {
// Volume of `a` that lies outside `b`. ~0 means a is (almost) contained in b.
double VolumeOutside(const Manifold& a, const Manifold& b) {
  return (a - b).Volume();
}
}  // namespace

// A purely convex shape has no concavities for any ball to fail to reach, so
// closing must leave it unchanged, and opening only changes it once the ball
// is larger than the feature (radius > shape radius).
TEST(Morphology, CloseConvexUnchanged) {
  const Manifold sphere = Manifold::Sphere(1.0, 32);
  const double v0 = sphere.Volume();
  for (bool aniso : {false, true}) {
    const Manifold closed = sphere.MorphologicalClose(0.5, 0, 50, aniso);
    EXPECT_EQ(closed.Status(), Manifold::Error::NoError);
    EXPECT_FALSE(closed.IsEmpty());
    EXPECT_NEAR(closed.Volume(), v0, 1e-6 * v0);
  }
}

// Opening with a ball smaller than the sphere reaches everywhere -> no change.
TEST(Morphology, OpenSmallRadiusUnchanged) {
  const Manifold sphere = Manifold::Sphere(1.0, 32);
  const double v0 = sphere.Volume();
  const Manifold opened = sphere.MorphologicalOpen(0.5, 0, 50, false);
  EXPECT_EQ(opened.Status(), Manifold::Error::NoError);
  EXPECT_NEAR(opened.Volume(), v0, 1e-6 * v0);
}

// A ball larger than the sphere fits nowhere inside it, so the morphological
// opening is (essentially) empty -- the convex surface flows all the way in.
TEST(Morphology, OpenLargerThanSphereVanishes) {
  const Manifold sphere = Manifold::Sphere(1.0, 64);
  const double v0 = sphere.Volume();
  const Manifold opened = sphere.MorphologicalOpen(4.0, 0.1, 200, false);
  EXPECT_EQ(opened.Status(), Manifold::Error::NoError);
  EXPECT_LT(opened.Volume(), 0.05 * v0);
}

// Opening shaves off a thin protrusion a ball cannot reach into, while the
// thick body the ball does fit inside is preserved.
TEST(Morphology, OpenRemovesProtrusion) {
  const Manifold body = Manifold::Cube({4, 4, 2}, true);
  const Manifold spike =
      Manifold::Cube({0.5, 0.5, 2.0}, true).Translate({0, 0, 2});
  const Manifold spiked = body + spike;
  const double vSpiked = spiked.Volume();

  const Manifold opened = spiked.MorphologicalOpen(0.5, 0.12, 60, true);
  EXPECT_EQ(opened.Status(), Manifold::Error::NoError);
  EXPECT_FALSE(opened.IsEmpty());
  EXPECT_LT(opened.Volume(), vSpiked);             // spike is shaved away
  EXPECT_GT(opened.Volume(), 0.5 * body.Volume());  // body survives
  EXPECT_LT(VolumeOutside(opened, spiked), 1e-2 * vSpiked);  // contained
}

// Closing a slot cut into a cube fills it in: volume grows but never exceeds
// the original solid cube, and the notched input stays contained in the result.
TEST(Morphology, CloseCubeNotch) {
  const Manifold cube = Manifold::Cube({2, 2, 2}, true);
  const double vCube = cube.Volume();
  // A narrow slot cut down into the top face.
  const Manifold cutter =
      Manifold::Cube({0.4, 2.2, 0.8}, true).Translate({0, 0, 0.8});
  const Manifold notched = cube - cutter;
  const double vNotched = notched.Volume();
  ASSERT_LT(vNotched, vCube);

  for (bool aniso : {false, true}) {
    const Manifold closed = notched.MorphologicalClose(0.6, 0.1, 60, aniso);
    EXPECT_EQ(closed.Status(), Manifold::Error::NoError);
    EXPECT_FALSE(closed.IsEmpty());
    // The slot gets (partially) filled: volume increases.
    EXPECT_GT(closed.Volume(), vNotched + 1e-3);
    // ...but closing never removes material beyond the solid cube.
    EXPECT_LT(closed.Volume(), vCube + 1e-3);
    // Closing contains the input.
    EXPECT_LT(VolumeOutside(notched, closed), 1e-2 * vNotched);
  }
}

// Closing the concave neck between two overlapping spheres adds material there.
TEST(Morphology, CloseTwoSphereNeck) {
  const Manifold a = Manifold::Sphere(1.0, 48).Translate({-0.6, 0, 0});
  const Manifold b = Manifold::Sphere(1.0, 48).Translate({0.6, 0, 0});
  const Manifold dumbbell = a + b;
  const double v0 = dumbbell.Volume();

  const Manifold closed = dumbbell.MorphologicalClose(1.5, 0.1, 60, true);
  EXPECT_EQ(closed.Status(), Manifold::Error::NoError);
  EXPECT_FALSE(closed.IsEmpty());
  EXPECT_GT(closed.Volume(), v0);                       // neck filled outward
  EXPECT_LT(VolumeOutside(dumbbell, closed), 1e-2 * v0);  // contains input
}

// Degenerate / no-op inputs are handled gracefully.
TEST(Morphology, DegenerateInputs) {
  const Manifold cube = Manifold::Cube({1, 1, 1}, true);
  // radius 0 is a no-op.
  const Manifold zero = cube.MorphologicalClose(0.0);
  EXPECT_NEAR(zero.Volume(), cube.Volume(), 1e-9);
  // Empty input stays empty.
  const Manifold empty = Manifold().MorphologicalClose(1.0);
  EXPECT_TRUE(empty.IsEmpty());
}
