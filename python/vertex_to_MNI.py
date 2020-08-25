#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 11:35:20 2020

@author: lau
"""

from scipy import linalg
import numpy as np
from functools import partial
from sklearn.neighbors import BallTree
from copy import deepcopy

def transform_coords(coords_in, tfm):
    """Affine coordinate transformations."""
    coords_out = np.dot(tfm[0:3, 0:3], coords_in.T) + tfm[0:3, 3]

    return coords_out

def stc_2_mgzvol(coords_in, fwd, mri_mgz):
    """Convert stc coordinates to match nilearn plotting on mgz MRI."""
    # convert to mm
    coords_in = deepcopy(coords_in)
    coords_in *= 1000

    # MEG headspace to RAS surface
    ras2meg = deepcopy(fwd['mri_head_t']['trans'])
    ras2meg[0:3, 3] *= 1000.

    coords_ras = transform_coords(coords_in, np.linalg.inv(ras2meg))

    # RAS surface to mgz voxel space
    vox2ras = mri_mgz.header.get_vox2ras_tkr()
    coords_mgz = transform_coords(coords_ras, np.linalg.inv(vox2ras))

    # mgz voxel space to world space
    mgz2mri = mri_mgz.header.get_affine()
    coords_out = transform_coords(coords_mgz, mgz2mri)

    return coords_out

class _DistanceQuery(object):
    """Wrapper for fast distance queries."""

    def __init__(self, xhs, method='BallTree', allow_kdtree=False):
        assert method in ('BallTree', 'cKDTree', 'cdist')

        # Fastest for our problems: balltree
        if method == 'BallTree':
            try:
                from sklearn.neighbors import BallTree
            except ImportError:
                # logger.info('Nearest-neighbor searches will be significantly '
                #             'faster if scikit-learn is installed.')
                method = 'cKDTree'
            else:
                self.query = partial(_safe_query, func=BallTree(xhs).query,
                                     reduce=True, return_distance=True)

        # Then cKDTree
        if method == 'cKDTree':
            try:
                from scipy.spatial import cKDTree
            except ImportError:
                method = 'cdist'
            else:
                self.query = cKDTree(xhs).query

        # KDTree is really only faster for huge (~100k) sets,
        # (e.g., with leafsize=2048), and it's slower for small (~5k)
        # sets. We can add it later if we think it will help.

        # Then the worst: cdist
        # if method == 'cdist':
        #     self.query = _CDist(xhs).query

        self.data = xhs

def apply_trans(trans, pts, move=True):
    """Apply a transform matrix to an array of points.

    Parameters
    ----------
    trans : array, shape = (4, 4) | instance of Transform
        Transform matrix.
    pts : array, shape = (3,) | (n, 3)
        Array with coordinates for one or n points.
    move : bool
        If True (default), apply translation.

    Returns
    -------
    transformed_pts : shape = (3,) | (n, 3)
        Transformed point(s).
    """
    if isinstance(trans, dict):
        trans = trans['trans']
    pts = np.asarray(pts)
    if pts.size == 0:
        return pts.copy()

    # apply rotation & scale
    out_pts = np.dot(pts, trans[:3, :3].T)
    # apply translation
    if move:
        out_pts += trans[:3, 3]

    return out_pts

def _safe_query(rr, func, reduce=False, **kwargs):
    if len(rr) == 0:
        return np.array([]), np.array([], int)
    out = func(rr)
    out = [out[0][:, 0], out[1][:, 0]] if reduce else out
    return out


def _cut_coords_to_ijk(cut_coords, img):
    ijk = apply_trans(linalg.inv(img.affine), cut_coords)
    ijk = np.clip(np.round(ijk).astype(int), 0, np.array(img.shape[:3]) - 1)
    return ijk

def _cut_coords_to_idx(cut_coords, img):
        """Convert voxel coordinates to index in stc.data."""
        ijk = _cut_coords_to_ijk(cut_coords, img)
        del cut_coords
        dist, loc_idx = dist_to_verts.query(ijk[np.newaxis])
        dist, loc_idx = dist[0], loc_idx[0]
        return loc_idx

pos = np.array([-0.016, -0.015, 0.061])
img = stc.as_volume(src, mri_resolution=False)
stc_ijk = np.array(
            np.unravel_index(stc.vertices[0], img.shape[:3], order='F')).T
ijk = _cut_coords_to_ijk(pos, img)
dist_to_verts = _DistanceQuery(stc_ijk, allow_kdtree=True)
dist, loc_idx = dist_to_verts.query(ijk[np.newaxis])
dist, loc_idx = dist[0], loc_idx[0]
ijk = stc_ijk[loc_idx]

vertex = stc.vertices[0][loc_idx]
pos_vertex = src[0]['rr'][vertex, :]





