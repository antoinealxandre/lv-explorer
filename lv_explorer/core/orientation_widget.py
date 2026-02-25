"""
Widget d'orientation 3D avec buste humain stylisé.
Utilise vtkOrientationMarkerWidget pour un affichage fiable en bas à droite,
synchronisé automatiquement avec la caméra, comme dans les logiciels professionnels
(ParaView, 3D Slicer, etc.).
"""

import numpy as np
import vtk


def _make_actor_from_source(source, color, position=(0, 0, 0),
                            scale=(1, 1, 1), ambient=0.3, diffuse=0.7,
                            specular=0.1):
    """Helper : crée un vtkActor depuis un vtkAlgorithm."""
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(source.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(*color)
    actor.GetProperty().SetAmbient(ambient)
    actor.GetProperty().SetDiffuse(diffuse)
    actor.GetProperty().SetSpecular(specular)
    actor.SetPosition(*position)
    actor.SetScale(*scale)
    return actor


def create_human_bust_prop():
    """
    Crée un buste humain en chargeant un fichier .obj depuis assets/human_bust.obj
    (compatible avec vtkOrientationMarkerWidget).

    Si le fichier n'existe pas, retourne un buste approximatif de secours.

    Orientation anatomique :
      X+ = droite patient (Right)
      Y+ = supérieur (Superior)
      Z+ = antérieur (Anterior)

    Returns
    -------
    vtk.vtkPropAssembly
    """
    import os
    from pathlib import Path
    from ..utils import resource_path

    assembly = vtk.vtkPropAssembly()

    # Chemin vers human_bust.obj — compatible développement ET bundle PyInstaller
    obj_path = Path(resource_path(os.path.join("assets", "human_bust.obj")))

    # Essayer de charger le fichier .obj
    if obj_path.exists():
        try:
            reader = vtk.vtkOBJReader()
            reader.SetFileName(str(obj_path))
            reader.Update()
            
            # Appliquer rotations pour aligner avec le système Y+=Superior, Z+=Anterior
            # Il faut redresser le buste ET le faire face à nous
            transform = vtk.vtkTransform()
            transform.RotateX(90)    # Redresser : Z → Y (haut)
            transform.RotateZ(180)   # Faire face : retourner de 180°
            
            transform_filter = vtk.vtkTransformPolyDataFilter()
            transform_filter.SetInputConnection(reader.GetOutputPort())
            transform_filter.SetTransform(transform)
            transform_filter.Update()
            
            # Créer un actor depuis le mesh transformé
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(transform_filter.GetOutputPort())
            
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetColor(0.85, 0.75, 0.65)  # Couleur peau
            actor.GetProperty().SetAmbient(0.3)
            actor.GetProperty().SetDiffuse(0.7)
            actor.GetProperty().SetSpecular(0.2)
            
            # Normaliser la taille (auto-scale pour que le buste soit dans [-1, 1])
            bounds = transform_filter.GetOutput().GetBounds()
            max_dim = max(bounds[1] - bounds[0], 
                         bounds[3] - bounds[2], 
                         bounds[5] - bounds[4])
            if max_dim > 0:
                scale = 1.6 / max_dim  # Ajuster la taille
                actor.SetScale(scale, scale, scale)
            
            # Centrer le modèle
            center_x = (bounds[0] + bounds[1]) / 2.0
            center_y = (bounds[2] + bounds[3]) / 2.0
            center_z = (bounds[4] + bounds[5]) / 2.0
            actor.SetPosition(-center_x * actor.GetScale()[0],
                            -center_y * actor.GetScale()[1],
                            -center_z * actor.GetScale()[2])
            
            assembly.AddPart(actor)
            print(f"✓ Loaded human bust from {obj_path}")
            
        except Exception as e:
            print(f"⚠ Error loading {obj_path}: {e}")
            print("  → Using fallback geometric bust")
            _create_fallback_bust(assembly)
    else:
        print(f"⚠ File not found: {obj_path}")
        print("  → Using fallback geometric bust")
        _create_fallback_bust(assembly)
    
    # Toujours ajouter les labels d'orientation
    _add_orientation_labels(assembly)
    
    return assembly


def _create_fallback_bust(assembly):
    """Crée un buste géométrique approximatif de secours.
    
    Le buste est construit dans l'orientation standard (Z+ = haut),
    puis tourné pour s'aligner avec le système où Y+ = Superior et Z+ = Anterior.
    Rotations : 90° autour de X (redresser) puis 180° autour de Z (faire face).
    """
    SKIN = (0.85, 0.75, 0.65)
    
    # Transformation globale
    transform = vtk.vtkTransform()
    transform.RotateX(90)     # Z → Y (redresser le buste)
    transform.RotateZ(180)    # Retourner pour faire face
    
    # Torse (ellipsoïde allongé verticalement en Z)
    torso = vtk.vtkParametricEllipsoid()
    torso.SetXRadius(0.35)
    torso.SetYRadius(0.20)
    torso.SetZRadius(0.45)
    torso_src = vtk.vtkParametricFunctionSource()
    torso_src.SetParametricFunction(torso)
    torso_src.SetUResolution(24)
    torso_src.SetVResolution(24)
    
    torso_filter = vtk.vtkTransformPolyDataFilter()
    torso_filter.SetInputConnection(torso_src.GetOutputPort())
    torso_filter.SetTransform(transform)
    torso_filter.Update()
    assembly.AddPart(_make_actor_from_source(torso_filter, SKIN, position=(0, -0.1, 0)))
    
    # Cou (cylindre vertical le long de Z)
    neck = vtk.vtkCylinderSource()
    neck.SetRadius(0.10)
    neck.SetHeight(0.15)
    neck.SetResolution(16)
    
    neck_filter = vtk.vtkTransformPolyDataFilter()
    neck_filter.SetInputConnection(neck.GetOutputPort())
    neck_filter.SetTransform(transform)
    neck_filter.Update()
    assembly.AddPart(_make_actor_from_source(neck_filter, SKIN, position=(0, 0.40, 0)))
    
    # Tête (sphère en haut)
    head = vtk.vtkSphereSource()
    head.SetRadius(0.18)
    head.SetThetaResolution(20)
    head.SetPhiResolution(20)
    
    head_filter = vtk.vtkTransformPolyDataFilter()
    head_filter.SetInputConnection(head.GetOutputPort())
    head_filter.SetTransform(transform)
    head_filter.Update()
    assembly.AddPart(_make_actor_from_source(head_filter, SKIN, position=(0, 0.58, 0)))
    
    # Épaules
    for side in (-1, 1):
        sh = vtk.vtkSphereSource()
        sh.SetRadius(0.12)
        sh.SetThetaResolution(16)
        sh.SetPhiResolution(16)
        
        # Échelle puis rotation
        t_scale = vtk.vtkTransform()
        t_scale.Scale(1.0, 0.8, 0.7)
        t_scale.Concatenate(transform)
        
        sf = vtk.vtkTransformPolyDataFilter()
        sf.SetInputConnection(sh.GetOutputPort())
        sf.SetTransform(t_scale)
        sf.Update()
        assembly.AddPart(_make_actor_from_source(sf, SKIN, position=(side * 0.35, 0.25, 0)))
    
    # Nez (cône pointant vers Y+ avant rotations, pour pointer vers Z+ après)
    nose = vtk.vtkConeSource()
    nose.SetRadius(0.04)
    nose.SetHeight(0.08)
    nose.SetResolution(12)
    nose.SetDirection(0, 1, 0)  # Pointe vers Y+ (deviendra Z+ après rotations)
    
    nose_filter = vtk.vtkTransformPolyDataFilter()
    nose_filter.SetInputConnection(nose.GetOutputPort())
    nose_filter.SetTransform(transform)
    nose_filter.Update()
    assembly.AddPart(_make_actor_from_source(nose_filter, (0.82, 0.68, 0.58), position=(0, 0.58, 0.19)))
    
    # Marqueur antérieur (devant)
    dot = vtk.vtkSphereSource()
    dot.SetRadius(0.04)
    dot.SetThetaResolution(10)
    dot.SetPhiResolution(10)
    
    dot_filter = vtk.vtkTransformPolyDataFilter()
    dot_filter.SetInputConnection(dot.GetOutputPort())
    dot_filter.SetTransform(transform)
    dot_filter.Update()
    assembly.AddPart(_make_actor_from_source(dot_filter, (0.8, 0.2, 0.2), position=(0, 0.05, 0.21), ambient=0.6))


def _add_orientation_labels(assembly):
    """Ajoute les labels R/L/A/P/S/I à l'assembly."""
    labels_def = [
        ('L', ( 0.70,  0,     0   ), (0.2, 0.4, 1.0)),  # Left (patient right = screen left)
        ('R', (-0.70,  0,     0   ), (0.2, 0.4, 1.0)),  # Right (patient left = screen right)
        ('S', ( 0,     0.85,  0   ), (1.0, 0.3, 0.1)),  # Superior (Y+)
        ('I', ( 0,    -0.68,  0   ), (1.0, 0.3, 0.1)),  # Inferior (Y-)
        ('A', ( 0,     0,     0.50), (0.1, 0.7, 0.1)),  # Anterior (Z+)
        ('P', ( 0,     0,    -0.50), (0.1, 0.7, 0.1)),  # Posterior (Z-)
    ]
    for text, pos, color in labels_def:
        vt = vtk.vtkVectorText()
        vt.SetText(text)
        vt.Update()
        bounds = vt.GetOutput().GetBounds()
        cx = -(bounds[0] + bounds[1]) / 2.0
        cy = -(bounds[2] + bounds[3]) / 2.0
        tc = vtk.vtkTransform()
        tc.Translate(cx, cy, 0)
        tf = vtk.vtkTransformPolyDataFilter()
        tf.SetInputConnection(vt.GetOutputPort())
        tf.SetTransform(tc)
        tf.Update()
        a = _make_actor_from_source(tf, color, position=pos,
                                    scale=(0.12, 0.12, 0.12),
                                    ambient=1.0, diffuse=0.0, specular=0.0)
        assembly.AddPart(a)


class HumanBustOrientationWidget:
    """
    Widget d'orientation avec buste humain, bas-droite de la fenêtre.
    Utilise vtkOrientationMarkerWidget — fiable, pas besoin de gestion
    manuelle de renderer/layer/camera.

    Usage
    -----
        widget = HumanBustOrientationWidget()
        widget.setup(interactor)
    """

    def __init__(self):
        self._widget = None

    def setup(self, interactor):
        """
        Initialise le widget sur le vtkRenderWindowInteractor donné.

        Parameters
        ----------
        interactor : vtkRenderWindowInteractor
            L'interactor de la fenêtre de rendu (plotter.iren.interactor
            pour pyvistaqt).
        """
        if interactor is None:
            return

        bust = create_human_bust_prop()

        self._widget = vtk.vtkOrientationMarkerWidget()
        self._widget.SetOrientationMarker(bust)
        self._widget.SetInteractor(interactor)
        self._widget.SetViewport(0.0, 0.0, 0.22, 0.22)   # bas-droite
        self._widget.EnabledOn()
        self._widget.InteractiveOff()   # pas de drag interactif

    def cleanup(self):
        """Supprime le widget proprement."""
        if self._widget is not None:
            try:
                self._widget.EnabledOff()
            except Exception:
                pass
            self._widget = None

    @property
    def enabled(self):
        return self._widget is not None
