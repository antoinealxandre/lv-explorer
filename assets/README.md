# Assets Directory

This folder contains 3D models and other assets for the LV Explorer application.

## Human Bust Model

Place your 3D human bust model here as `human_bust.obj`.

**Expected file:** `human_bust.obj`

### Requirements:

- Format: Wavefront OBJ (.obj)
- Orientation:
  - X+ = Right (patient's right side)
  - Y+ = Anterior (front)
  - Z+ = Superior (up)
- The model will be automatically centered and scaled to fit the orientation widget

### Fallback Behavior:

If `human_bust.obj` is not found, the application will use a simplified geometric bust made of basic shapes (spheres, cylinders, ellipsoids).

## Adding the Model

1. Export or download your human bust 3D model in OBJ format
2. Name it `human_bust.obj`
3. Place it in this `assets/` directory
4. Restart the application to see the new model in the orientation widget (bottom-right corner)
