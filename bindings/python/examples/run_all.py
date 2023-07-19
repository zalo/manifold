"""
 Copyright 2022 The Manifold Authors.

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

      https://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 """

import pathlib
import sys
import importlib
import trimesh
import numpy as np
from time import time

if __name__ == "__main__":
    current_file = pathlib.Path(__file__)
    current_dir = current_file.parent
    files = [f.parts[-1][:-3]
             for f in current_dir.glob('*.py') if f != current_file]

    export_models = len(sys.argv) == 2 and sys.argv[-1] == '-e'

    for f in files:
        module = importlib.import_module(f)
        t0 = time()
        model = module.run()
        mesh = model.to_meshgl()
        if export_models:
            vertices = np.reshape(mesh.vert_properties, (-1, mesh.num_prop))
            if mesh.num_prop > 3:
                vertices = np.hsplit(vertices, 3)
            meshOut = trimesh.Trimesh(
                vertices=vertices, faces=np.reshape(mesh.tri_verts, (-1, 3)))
            trimesh.exchange.export.export_mesh(meshOut, f'{f}.glb', 'glb')
            print(f'Exported model to {f}.glb')
        t1 = time()
        print(f'Took {(t1-t0)*1000:.1f}ms for {f}')
