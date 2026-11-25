import numpy as np
import cv2
import random
from scipy.interpolate import splprep, splev
from skimage.filters import sobel, gaussian
from skimage.draw import line as draw_line_pixels
from skimage.morphology import closing, opening, disk, dilation, erosion, binary_closing, binary_opening
from scipy.ndimage import distance_transform_edt
from typing import List, Tuple, Dict
import matplotlib
# matplotlib.use('Agg')  # Force non-interactive backend before importing pyplot
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from PIL import Image
import os
import time
from tqdm import tqdm

################################################################################
#                           UTILITY FUNCTIONS                                  #
################################################################################

def compute_max_buffer_distance(width: int, height: int, num_neurons: int,
                                mean_soma_radius: float, std_soma_radius: float,
                                edge_margin: int) -> int:
    """Greedy search for a square packing buffer so all somas fit on the canvas."""
    max_soma_radius = mean_soma_radius + 2 * std_soma_radius
    eff_w           = width  - 2 * edge_margin
    eff_h           = height - 2 * edge_margin
    low, high       = 0, max(0, min(eff_w, eff_h) - 2 * max_soma_radius)
    best            = 0
    while low <= high:
        mid = (low + high) // 2
        cell = 2 * max_soma_radius + mid
        cols = int(eff_w // cell)
        rows = int(eff_h // cell)
        if cols > 0 and rows > 0 and cols * rows >= num_neurons:
            best, low = mid, mid + 1
        else:
            high = mid - 1
    return int(best)

_DISC_CACHE = {}
def _disc_offsets(radius: int):
    """Integer offsets inside a disc (origin excluded) as float arrays (dc, dr, norm), cached."""
    radius = int(radius)
    if radius not in _DISC_CACHE:
        dr, dc = np.mgrid[-radius:radius + 1, -radius:radius + 1]
        keep = (dc * dc + dr * dr <= radius * radius) & ((dc != 0) | (dr != 0))
        dc = dc[keep].astype(float); dr = dr[keep].astype(float)
        _DISC_CACHE[radius] = (dc, dr, np.sqrt(dc * dc + dr * dr))
    return _DISC_CACHE[radius]

################################################################################
#                               SOMA CLASS                                     #
################################################################################

class Soma:
    """Blob‑shaped soma, drawn once and then treated as an immutable disc."""
    def __init__(self, position: Tuple[float, float], mean_r: float, std_r: float):
        self.position = position
        self.radius   = max(np.random.normal(mean_r, std_r), 1.0)
        self._generate()

    # ------------------------------------------------------------------ helpers
    def _generate(self):
        theta     = np.linspace(0, 2 * np.pi, 100)
        mask_size = max(int(self.radius * 4), 5)
        mask      = np.zeros((mask_size, mask_size), np.uint8)
        c         = mask_size // 2

        def add_blob(cx, cy, r, ratio=1.0, rot=0.0):
            pts = np.column_stack([
                cx + r * ( np.cos(theta) * np.cos(rot) - ratio * np.sin(theta) * np.sin(rot)),
                cy + r * ( np.sin(theta) * np.cos(rot) + ratio * np.cos(theta) * np.sin(rot))
            ]).astype(np.int32).reshape(-1, 1, 2)
            cv2.fillPoly(mask, [pts], 1)

        add_blob(c, c, 0.8 * self.radius)
        for _ in range(np.random.randint(7, 13)):
            a     = np.random.uniform(0, 2 * np.pi)
            d     = np.random.uniform(0, 0.7 * self.radius)
            cx, cy= c + d * np.cos(a), c + d * np.sin(a)
            r     = np.random.uniform(0.4, 0.8) * self.radius
            add_blob(cx, cy, r, ratio=np.random.uniform(0.8, 1.2), rot=np.random.uniform(0, 2*np.pi))

        if np.any(mask):
            mask = cv2.GaussianBlur(mask, (5,5), 0)
        cnt, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        if cnt and len(cnt[0]) > 0:
            cnt_squeezed = cnt[0].squeeze()
            if cnt_squeezed.ndim == 1 and cnt_squeezed.shape[0] == 2:
                 self.x = np.array([cnt_squeezed[0] - c + self.position[0]])
                 self.y = np.array([cnt_squeezed[1] - c + self.position[1]])
            elif cnt_squeezed.ndim == 2:
                 self.x = cnt_squeezed[:,0] - c + self.position[0]
                 self.y = cnt_squeezed[:,1] - c + self.position[1]
            else:
                 th    = np.linspace(0, 2*np.pi, 32)
                 self.x= self.position[0] + self.radius * np.cos(th)
                 self.y= self.position[1] + self.radius * np.sin(th)
        else:
            th    = np.linspace(0, 2*np.pi, 32)
            self.x= self.position[0] + self.radius * np.cos(th)
            self.y= self.position[1] + self.radius * np.sin(th)
        self.local_mask = (mask>0).astype(np.uint8)

    # ----------------------------------------------------------- public helpers
    def binary_masks(self, canvas_shape: Tuple[int,int], label:int=1):
        pts = np.column_stack([self.x, self.y]).astype(np.int32)
        filled   = np.zeros(canvas_shape, np.uint16 if label>1 else np.uint8)
        skeleton = np.zeros_like(filled)
        cv2.fillPoly(filled,   [pts], int(label))
        cv2.polylines(skeleton,[pts], True,  int(label), 1)
        return filled, skeleton

################################################################################
#                             DENDRITE CLASS                                   #
################################################################################

class Dendrite:
    """Holds global parameters and branch‑length schedule."""
    def __init__(self, soma:Soma, depth:int, D:float, mean_br:int, branch_angle:float,
                 weave_type:str='Gauss', randomness:float=0.1, curviness:str='Gauss',
                 curviness_mag:float=1.0, n_primary:int=4, total_len:float=500,
                 init_th:float=10, deterministic=False):
        self.soma            = soma
        self.depth           = depth
        self.D               = D if D > 0 else 1.0
        self.mean_branches   = mean_br if mean_br > 0 else 1
        self.branch_angle    = branch_angle
        self.weave_type      = weave_type
        self.randomness      = randomness
        self.curviness       = curviness
        self.curviness_mag   = curviness_mag
        self.n_primary       = n_primary
        self.total_length    = total_len if total_len > 0 else 1.0
        self.initial_thick   = init_th
        self.det             = deterministic
        # schedule of straight‑line lengths per depth level
        if self.depth > 0:
            r_val = self.mean_branches ** (-1 / self.D) if self.mean_branches > 0 and self.D > 0 else 0
            sum_r_i = sum(r_val**i for i in range(self.depth))
            if sum_r_i == 0: sum_r_i = 1e-6
            norm = self.total_length / sum_r_i
            self.branch_lengths = np.array([norm * r_val**i for i in range(self.depth)])
        else:
            self.branch_lengths = np.array([])

    def intra_branch_weave(self, x1: float, y1: float, x2: float, y2: float, length: float):
        """
        Creates a curved path between two points using splines, similar to bn_new.py.
        
        Args:
            x1, y1: Start point coordinates
            x2, y2: End point coordinates  
            length: Length of the segment
            
        Returns:
            tuple: (xs, ys) arrays of points along the curved path
        """
        if self.det:
            return np.array([x1, x2]), np.array([y1, y2])
        
        # Parameters for curve generation
        num_points = max(int(self.curviness_mag * 10), 2)
        num_control_points = np.random.randint(4, 7)

        # Create parameter values for control points
        t_values = np.linspace(0.2, 0.8, num_control_points - 2)
        t_values = np.concatenate(([0], t_values, [1]))

        # Calculate the base direction and perpendicular angle
        dx = x2 - x1
        dy = y2 - y1
        branch_angle_rad = np.arctan2(dy, dx)
        perp_angle_rad = branch_angle_rad + np.pi / 2

        # Create base points along the straight line
        base_x = x1 + t_values * dx
        base_y = y1 + t_values * dy

        # Add perpendicular perturbations to create curvature
        perturbation_magnitude = length * np.random.uniform(-0.05, 0.05, size=len(t_values))
        control_x = base_x + perturbation_magnitude * np.cos(perp_angle_rad)
        control_y = base_y + perturbation_magnitude * np.sin(perp_angle_rad)

        # Create a spline through the control points
        try:
            tck, u = splprep([control_x, control_y], s=0)
            u_fine = np.linspace(0, 1, num_points)
            xs, ys = splev(u_fine, tck)
        except:
            # Fallback to straight line if spline fails
            xs = np.linspace(x1, x2, num_points)
            ys = np.linspace(y1, y2, num_points)

        # Add fine-scale noise for wiggle
        perturbation_scale = length / 200
        if self.curviness == 'Gauss':
            xs += np.random.normal(0, perturbation_scale, num_points)
            ys += np.random.normal(0, perturbation_scale, num_points)
        elif self.curviness == 'Uniform':
            xs += np.random.uniform(-perturbation_scale, perturbation_scale, num_points)
            ys += np.random.uniform(-perturbation_scale, perturbation_scale, num_points)

        # Ensure endpoints are preserved
        xs[0], ys[0] = x1, y1
        xs[-1], ys[-1] = x2, y2

        return xs, ys

################################################################################
#                             NEURON CLASS                                     #
################################################################################

class Neuron:
    def __init__(self, position:Tuple[float,float], params:dict, network:'Network', nid:str):
        self.network   = network
        self.position  = position
        self.start_a   = position[2] if len(position)==3 else 0.0
        
        # Create soma first to get its actual radius
        soma = Soma(position[:2], params['mean_soma_radius'], params['std_soma_radius'])
        
        # Calculate initial thickness proportional to soma size
        if 'thickness_soma_ratio' in params:
            # Use proportional thickness based on soma radius
            calculated_thickness = soma.radius * params['thickness_soma_ratio']
            # Cap at maximum thickness for reasonable dendrite sizes
            calculated_thickness = min(calculated_thickness, 30.0)
        else:
            # Fallback to fixed thickness for backward compatibility
            calculated_thickness = params.get('initial_thickness', 10)
            # Also cap fixed thickness at maximum
            calculated_thickness = min(calculated_thickness, 30.0)
        
        self.dendrite  = Dendrite(soma,
                                  params['depth'], params.get('D', 2.0), params.get('mean_branches',3),
                                  params.get('branch_angle', np.pi/4), params.get('weave_type','Gauss'),
                                  params.get('randomness',0.1), params.get('curviness','Gauss'),
                                  params.get('curviness_magnitude',1.0), params.get('n_primary_dendrites',4),
                                  params.get('total_length',500), calculated_thickness,
                                  network.deterministic)
        self.nid_str   = nid
        self.id_int    = None
        self.pass_prob = params.get('pass_through_probability',0.0)
        self.fusion_prob = params.get('fusion_probability',0.0)
        self.fusion_distance_threshold = params.get('fusion_distance_threshold',10.0)
        self.is_growing= True
        self.branch_data : List[Dict]= []
        self.visible_total_length = float(self.dendrite.total_length)
        self.max_recorded_path_length = 0.0
        self.connections: List[Tuple[float, float, int]] = []  # Track connection points (x, y, step)
        self.fusion_points: List[Tuple[float, float, int, int]] = []  # Track fusion points (x, y, other_neuron_id, step)
        self.soma_connections: List[Tuple[float, float, int, int]] = []  # Track soma-dendrite connections (x, y, soma_id, step)
        
        # Spline curvature parameters
        self.spline_segment_length = params.get('spline_segment_length', 25.0)  # Length over which to apply splines
        self.use_splines = params.get('use_splines', True)  # Whether to use spline curvature
        
        # Trajectory planning parameters
        self.trajectory_planning = params.get('trajectory_planning', True)  # Enable trajectory-based growth
        self.trajectory_length = params.get('trajectory_length', 50.0)  # Distance ahead to plan trajectories
        self.trajectory_deviation = params.get('trajectory_deviation', 15.0)  # How much trajectories can curve
        self.trajectory_replan_interval = params.get('trajectory_replan_interval', 20.0)  # When to replan trajectories
        self.trajectory_smoothness = params.get('trajectory_smoothness', 0.1)  # Spline smoothing parameter
        
        # Self-avoidance (raster-local, low-overhead)
        self.self_avoidance = params.get('self_avoidance', False)
        self.self_avoidance_radius = params.get('self_avoidance_radius', 4)
        self.self_avoidance_retries = params.get('self_avoidance_retries', 3)
        self.self_avoidance_parent_exempt_radius = params.get('self_avoidance_parent_exempt_radius', 4)
        self.self_avoidance_min_path_distance = params.get('self_avoidance_min_path_distance', 12.0)
        self.self_avoidance_cone_half_angle = params.get('self_avoidance_cone_half_angle', np.pi / 2)
        
        # Variable primary dendrites parameters
        self.variable_primary_dendrites = params.get('variable_primary_dendrites', True)  # Enable variable primary dendrites
        self.primary_dendrites_std = params.get('primary_dendrites_std', 1.0)  # Standard deviation for Gaussian sampling
        
        # Variable growth rate parameters
        base_growth_rate = params.get('growth_rate', 2.5)  # Base growth rate from parameters
        self.variable_growth_rate = params.get('variable_growth_rate', True)  # Enable variable growth rates
        growth_rate_std = params.get('growth_rate_std', 0.5)  # Standard deviation for growth rate sampling
        
        # Sample individual growth rate for this neuron
        if self.variable_growth_rate and not self.dendrite.det:
            # Sample from Gaussian with mean=base_growth_rate and configurable std
            sampled_growth_rate = np.random.normal(base_growth_rate, growth_rate_std)
            # Ensure positive growth rate with minimum of 0.1
            self.growth_rate = max(0.1, sampled_growth_rate)
        else:
            # Use fixed growth rate (deterministic or disabled variable rates)
            self.growth_rate = base_growth_rate
        
        # tip list for time integration
        self.active_tips: List[Dict] = []
        self._tip_id    = 0
        self.last_step_stats = {
            'self_avoidance_attempts': 0,
            'self_avoidance_rejections': 0,
            'tips_processed': 0,
            'growth_direction_select_s': 0.0,
            'growth_self_avoidance_s': 0.0,
            'growth_collision_scan_s': 0.0,
            'growth_draw_record_s': 0.0,
        }
        self._init_primary_tips()

    # -------------------------------------------------------- tip initialisation
    def _init_primary_tips(self):
        if self.dendrite.depth == 0:
            self.is_growing = False
            return

        xs, ys = self.dendrite.soma.x, self.dendrite.soma.y
        cx, cy = self.position[:2]
        
        # Sample number of primary dendrites from Gaussian distribution
        if self.variable_primary_dendrites and not self.dendrite.det:
            # Sample from Gaussian with mean=n_primary and configurable std
            sampled_n_primary = np.random.normal(self.dendrite.n_primary, self.primary_dendrites_std)
            # Ensure at least 1 dendrite and round to integer
            actual_n_primary = max(1, int(np.round(sampled_n_primary)))
        else:
            # Use fixed number (deterministic or disabled variable dendrites)
            actual_n_primary = self.dendrite.n_primary
        
        base_a = np.linspace(0, 2*np.pi, actual_n_primary, endpoint=False)
        base_a = (base_a + self.start_a) % (2*np.pi)
        if not self.dendrite.det:
            base_a += np.random.uniform(-np.pi/6, np.pi/6, size=base_a.shape)
        for a in base_a:
            # find contour point closest in angle
            angs = (np.arctan2(ys-cy, xs-cx)+2*np.pi)%(2*np.pi)
            idx  = np.argmin(np.abs((angs-a+np.pi)%(2*np.pi)-np.pi))
            pt   = np.array([xs[idx], ys[idx]], float)
            if not hasattr(self, '_primary_tip_templates'):
                self._primary_tip_templates = {}
            self._primary_tip_templates[self._tip_id] = (pt.copy(), float(a))
            self.active_tips.append(dict(id=self._tip_id, pt=pt, angle=a, depth=0, traced_len=0.0, 
                                        parent_tip_id=None,
                                        parent_branch_path_dist=0.0,
                                        encountered_dendrites=set(),
                                        pass_through_dendrites=set(),
                                        accumulated_points=[], accumulated_distances=0.0,
                                        accumulated_start_traced_len=0.0,
                                        # Trajectory planning data
                                        trajectory_points=None,  # List of target points along planned trajectory
                                        trajectory_index=0,      # Current position in trajectory
                                        last_replan_distance=0.0  # Distance since last trajectory replan
                                        ))
            self._tip_id+=1
        
        if not self.active_tips:
             self.is_growing = False

    # ------------------------------------------------------ growth helpers ----
    def _segment_thickness(self, depth:int, traced:float, seg_len_of_current_depth_level:float):
        # Calculate total distance from soma to the START of the current depth's segment
        dist_to_start_of_current_depth_level = 0
        if depth > 0:
            dist_to_start_of_current_depth_level = np.sum(self.dendrite.branch_lengths[:depth])

        prop_s = 1 - (dist_to_start_of_current_depth_level + traced) / self.dendrite.total_length
        prop_e = 1 - (dist_to_start_of_current_depth_level + traced + seg_len_of_current_depth_level) / self.dendrite.total_length
        
        prop_s = np.clip(prop_s, 0,1)
        prop_e = np.clip(prop_e, 0,1)

        safe_D = self.dendrite.D if self.dendrite.D > 0 else 1.0

        t0 = max(1, self.dendrite.initial_thick * prop_s**(1/safe_D))
        t1 = max(1, self.dendrite.initial_thick * prop_e**(1/safe_D))
        return t0, t1

    def _update_branch_length_schedule(self, new_total_length: float):
        """Update the dendrite length schedule without resetting growth state."""
        fixed = getattr(self, 'fixed_schedule_length', None)
        self.visible_total_length = float(max(1.0, float(new_total_length)))
        if fixed is not None:
            # Branching schedule is pinned to the neuron's intrinsic arbor size; the
            # visible/target length only decides how far along that arbor is realised.
            new_total_length = fixed
        self.dendrite.total_length = max(1.0, float(new_total_length))
        if self.dendrite.depth <= 0:
            self.dendrite.branch_lengths = np.array([])
            return
        safe_branches = max(float(self.dendrite.mean_branches), 1.0)
        safe_D = self.dendrite.D if self.dendrite.D > 0 else 1.0
        r_val = safe_branches ** (-1.0 / safe_D)
        terms = np.array([r_val ** i for i in range(self.dendrite.depth)], dtype=float)
        denom = float(np.sum(terms)) if np.sum(terms) > 0 else 1.0
        self.dendrite.branch_lengths = (self.dendrite.total_length / denom) * terms

    def set_visible_total_length(self, new_total_length: float):
        """
        Set the visible/pruned dendrite length threshold for mask replay.
        Does not delete recorded segment history.
        """
        self._update_branch_length_schedule(new_total_length)
        self.visible_total_length = float(max(0.0, new_total_length))

    def _draw_segment(self, p0:np.ndarray, p1:np.ndarray, t:float):
        if getattr(self.network, 'incremental_branch_mode', False):
            # The segment was already rasterised into the per-neuron count layers by
            # _record_branch_segment and composed into the masks (max-id ownership).
            # Drawing again here would overwrite ownership by draw order (and with a
            # different width rounding), leaving masks inconsistent with the counts.
            return
        p0i = tuple(np.round(p0).astype(int))
        p1i = tuple(np.round(p1).astype(int))
        
        # Check if both points are within boundaries - don't draw if either is outside
        h_net, w_net = self.network.filled_mask.shape
        if (0 <= p0i[0] < w_net and 0 <= p0i[1] < h_net and 
            0 <= p1i[0] < w_net and 0 <= p1i[1] < h_net):
            cv2.line(self.network.filled_mask,   p0i, p1i, self.id_int, thickness=int(max(1,t)))
            cv2.line(self.network.skeleton_mask, p0i, p1i, self.id_int, thickness=1)

    def _generate_trajectory(self, start_point: np.ndarray, start_angle: float, target_length: float) -> List[np.ndarray]:
        """
        Generate a smooth spline-based trajectory for dendrite growth.
        
        Args:
            start_point: Current tip position [x, y]
            start_angle: Current growth direction in radians
            target_length: Distance ahead to plan trajectory
            
        Returns:
            List of points along the planned trajectory
        """
        if not self.trajectory_planning:
            # Fallback to straight line
            end_point = start_point + target_length * np.array([np.cos(start_angle), np.sin(start_angle)])
            return [start_point.copy(), end_point]
        
        # Create control points for spline trajectory
        num_control_points = max(3, int(target_length / 15))  # More control points for longer trajectories
        
        # Base direction from current angle
        base_direction = np.array([np.cos(start_angle), np.sin(start_angle)])
        
        # Generate control points with increasing distance and curvature
        control_points = [start_point.copy()]
        
        for i in range(1, num_control_points):
            progress = i / (num_control_points - 1)  # 0 to 1
            
            # Base position along straight line
            base_pos = start_point + progress * target_length * base_direction
            
            # Add perpendicular deviation for curvature
            perp_direction = np.array([-base_direction[1], base_direction[0]])  # Perpendicular vector
            
            # Gradually increase deviation potential with distance
            max_deviation = self.trajectory_deviation * progress * np.random.uniform(0.3, 1.0)
            
            # Random deviation from straight path
            deviation = np.random.uniform(-max_deviation, max_deviation)
            
            # Add some forward/backward variation too
            forward_variation = np.random.uniform(-target_length * 0.1, target_length * 0.1)
            
            curved_pos = base_pos + deviation * perp_direction + forward_variation * base_direction
            
            # Ensure we don't go outside network boundaries
            curved_pos[0] = np.clip(curved_pos[0], 10, self.network.width - 10)
            curved_pos[1] = np.clip(curved_pos[1], 10, self.network.height - 10)
            
            control_points.append(curved_pos)
        
        # Create smooth spline through control points
        try:
            if len(control_points) >= 3:
                # Prepare data for spline fitting
                control_points_array = np.array(control_points)
                x_coords = control_points_array[:, 0]
                y_coords = control_points_array[:, 1]
                
                # Fit spline
                tck, u = splprep([x_coords, y_coords], s=self.trajectory_smoothness * target_length, k=min(3, len(control_points)-1))
                
                # Generate fine points along spline
                num_trajectory_points = max(10, int(target_length / 3))  # Point every ~3 pixels
                u_new = np.linspace(0, 1, num_trajectory_points)
                trajectory_x, trajectory_y = splev(u_new, tck)
                
                trajectory_points = [np.array([x, y]) for x, y in zip(trajectory_x, trajectory_y)]
                return trajectory_points
            else:
                # Fallback for too few control points
                end_point = start_point + target_length * base_direction
                return [start_point.copy(), end_point]
                
        except Exception:
            # Fallback if spline fitting fails
            end_point = start_point + target_length * base_direction
            return [start_point.copy(), end_point]
    
    def _get_next_trajectory_target(self, tip_dict: Dict) -> Tuple[np.ndarray, float]:
        """
        Get the next target point along the planned trajectory.
        
        Args:
            tip_dict: Tip data dictionary
            
        Returns:
            Tuple of (target_point, target_angle) for next growth step
        """
        current_pos = tip_dict['pt']
        
        # Check if we need to generate or regenerate trajectory
        needs_new_trajectory = (
            tip_dict['trajectory_points'] is None or 
            tip_dict['trajectory_index'] >= len(tip_dict['trajectory_points']) - 1 or
            tip_dict['last_replan_distance'] >= self.trajectory_replan_interval
        )
        
        if needs_new_trajectory:
            # Generate new trajectory
            current_angle = tip_dict['angle']
            tip_dict['trajectory_points'] = self._generate_trajectory(
                current_pos, current_angle, self.trajectory_length
            )
            tip_dict['trajectory_index'] = 0
            tip_dict['last_replan_distance'] = 0.0
        
        # Find next point on trajectory that's ahead of current position
        trajectory_points = tip_dict['trajectory_points']
        target_index = tip_dict['trajectory_index']
        
        # Look for next point that's far enough ahead
        min_step_distance = 2.0  # Minimum distance to next target
        
        for i in range(target_index, len(trajectory_points)):
            candidate_target = trajectory_points[i]
            distance = np.linalg.norm(candidate_target - current_pos)
            
            if distance >= min_step_distance:
                tip_dict['trajectory_index'] = i
                
                # Calculate angle to target
                direction = candidate_target - current_pos
                target_angle = np.arctan2(direction[1], direction[0])
                
                return candidate_target, target_angle
        
        # If we're at the end of trajectory, use current angle for straight growth
        current_angle = tip_dict['angle']
        target_pos = current_pos + 5.0 * np.array([np.cos(current_angle), np.sin(current_angle)])
        return target_pos, current_angle

    def _check_pixel_collision(self, p:Tuple[int,int], thickness:float=1.0):
        """
        Two-stage collision detection:
        1. Check if we're in another object's area (filled mask)
        2. Check if our centerline touches their centerline (skeleton mask)
        
        Returns: (collision_type, colliding_id, skeleton_contact) or None
        skeleton_contact: True if skeletons are touching, False if only filled areas overlap
        """
        c_coord, r_coord = p[0],p[1]
        if not (0<=r_coord<self.network.height and 0<=c_coord<self.network.width):
            return ('boundary', 0, True)
        
        # Check a small radius around the point to account for thickness
        check_radius = max(1, int(thickness // 2))
        
        # Track what we find in the area
        found_filled_collision = None
        found_skeleton_collision = None
        
        for dr in range(-check_radius, check_radius + 1):
            for dc in range(-check_radius, check_radius + 1):
                check_r = r_coord + dr
                check_c = c_coord + dc
                
                if not (0<=check_r<self.network.height and 0<=check_c<self.network.width):
                    continue
                
                # Check soma collision (always blocks) - prioritize soma
                sid = self.network.somas_mask[check_r,check_c]
                if sid != 0 and sid!=self.id_int:
                    # For somas, check if we're touching the soma boundary (treated as skeleton contact)
                    found_skeleton_collision = ('soma', sid)
                    break
                
                # Check dendrite skeleton collision (centerline to centerline)
                skeleton_id = self.network.skeleton_mask[check_r,check_c]
                if skeleton_id != 0 and skeleton_id != self.id_int:
                    found_skeleton_collision = ('dendrite', skeleton_id)
                
                # Check dendrite filled area collision (we're in their territory but not touching centerlines yet)
                filled_id = self.network.filled_mask[check_r,check_c]
                if filled_id != 0 and filled_id != self.id_int and found_filled_collision is None:
                    found_filled_collision = ('dendrite', filled_id)
        
        # Return skeleton contact if found (this triggers actual connection/fusion)
        if found_skeleton_collision is not None:
            collision_type, colliding_id = found_skeleton_collision
            return (collision_type, colliding_id, True)
        
        # Return filled contact if found (this allows continued growth until skeleton contact)
        if found_filled_collision is not None:
            collision_type, colliding_id = found_filled_collision
            return (collision_type, colliding_id, False)
            
        return None

    def _check_soma_boundary_contact(self, p: Tuple[int, int], soma_id: int) -> bool:
        """
        Check if point is touching the boundary of the specified soma.
        For dendrite-soma connections, we want contact with soma boundary.
        """
        c_coord, r_coord = p[0], p[1]
        
        # Check a small area around the point
        for dr in range(-2, 3):
            for dc in range(-2, 3):
                check_r = r_coord + dr
                check_c = c_coord + dc
                
                if not (0<=check_r<self.network.height and 0<=check_c<self.network.width):
                    continue
                    
                # If this pixel belongs to the soma, we're touching it
                if self.network.somas_mask[check_r, check_c] == soma_id:
                    return True
                    
        return False

    def _self_avoidance_conflict(self, p:Tuple[int,int], tip_dict:Dict, thickness:float,
                                 current_dist_from_soma: float, growth_dir: np.ndarray) -> bool:
        """
        Raster-local self-avoidance check against this neuron's existing arbor.
        Uses a small search radius and exempts the immediate parent/start neighborhood.
        """
        if (not self.self_avoidance) or current_dist_from_soma < self.self_avoidance_min_path_distance:
            return False

        c_coord, r_coord = p[0], p[1]
        if not (0 <= r_coord < self.network.height and 0 <= c_coord < self.network.width):
            return True

        check_radius = max(int(self.self_avoidance_radius), int(max(1, thickness // 2)))
        exempt_radius = max(0, int(self.self_avoidance_parent_exempt_radius))
        tip_pt = np.array(tip_dict.get('pt', [c_coord, r_coord]), dtype=float)
        g = np.array(growth_dir, dtype=float)
        g_norm = float(np.linalg.norm(g))
        if g_norm <= 1e-12:
            return False
        g = g / g_norm
        cone_half_angle = float(np.clip(self.self_avoidance_cone_half_angle, 1e-3, np.pi))
        cone_cos_threshold = float(np.cos(cone_half_angle))

        dc, dr, rel_norm = _disc_offsets(check_radius)   # disc offsets, origin excluded (cached)
        # Directional cone: only consider self pixels within the forward cone.
        cos_angle = (dc * g[0] + dr * g[1]) / rel_norm
        check_r = r_coord + dr.astype(np.int64)
        check_c = c_coord + dc.astype(np.int64)
        ok = ((cos_angle > cone_cos_threshold)
              & (check_r >= 0) & (check_r < self.network.height)
              & (check_c >= 0) & (check_c < self.network.width))
        if not ok.any():
            return False
        check_r = check_r[ok]; check_c = check_c[ok]
        hit = self.network.filled_mask[check_r, check_c] == self.id_int
        if not hit.any():
            return False
        if exempt_radius > 0:
            # Exempt immediate neighborhood around current tip origin to avoid blocking
            # continuation along the just-grown parent branch.
            d_tip_sq = (check_c[hit] - tip_pt[0])**2 + (check_r[hit] - tip_pt[1])**2
            return bool(np.any(d_tip_sq > exempt_radius * exempt_radius))
        return True

    def step(self, dt:float, current_step:int=0):
        if not self.active_tips:
            self.last_step_stats = {
                'self_avoidance_attempts': 0,
                'self_avoidance_rejections': 0,
                'tips_processed': 0,
                'growth_direction_select_s': 0.0,
                'growth_self_avoidance_s': 0.0,
                'growth_collision_scan_s': 0.0,
                'growth_draw_record_s': 0.0,
            }
            self.is_growing=False; return
        self.last_step_stats = {
            'self_avoidance_attempts': 0,
            'self_avoidance_rejections': 0,
            'tips_processed': 0,
            'growth_direction_select_s': 0.0,
            'growth_self_avoidance_s': 0.0,
            'growth_collision_scan_s': 0.0,
            'growth_draw_record_s': 0.0,
        }
        new_tips=[]
        segments_to_draw = []  # Store segments to draw after collision detection
        splined_segments_to_draw = []  # Store accumulated segments for spline drawing
        
        for tip_dict in list(self.active_tips):
            self.last_step_stats['tips_processed'] += 1
            depth = tip_dict['depth']
            if depth>=self.dendrite.depth:
                self.active_tips.remove(tip_dict); continue
            
            # Target length for the current segment at this depth
            seg_len_for_this_depth = self.dendrite.branch_lengths[depth]
            remain_in_segment  = seg_len_for_this_depth - tip_dict['traced_len']
            
            # Growth step for this iteration, capped by remaining length in segment
            base_growth_step = dt * self.growth_rate  # Use individual neuron's growth rate
            # Limit max growth step to prevent jumping over dendrites
            max_step_size = 2.0  # Maximum pixels per step to ensure collision detection
            current_growth_step = min(base_growth_step, remain_in_segment, max_step_size)
            
            if current_growth_step < 1e-3:
                if remain_in_segment < 1e-3:
                    # Draw accumulated splined segment before spawning children
                    if len(tip_dict['accumulated_points']) > 1:
                        avg_thickness = tip_dict.get('avg_accumulated_thickness', 1.0)
                        dist_to_start_of_depth_level = sum(self.dendrite.branch_lengths[:depth]) if depth > 0 else 0
                        chunk_start_abs = dist_to_start_of_depth_level + tip_dict.get('accumulated_start_traced_len', 0.0)
                        chunk_end_abs = chunk_start_abs + tip_dict.get('accumulated_distances', 0.0)
                        splined_segments_to_draw.append((
                            tip_dict['accumulated_points'].copy(),
                            avg_thickness,
                            float(chunk_start_abs),
                            float(chunk_end_abs),
                            dict(
                                tip_id=tip_dict['id'],
                                parent_tip_id=tip_dict.get('parent_tip_id'),
                                parent_branch_path_dist=float(tip_dict.get('parent_branch_path_dist', 0.0)),
                                depth=depth,
                                tip_traced_start=float(tip_dict.get('accumulated_start_traced_len', 0.0)),
                                tip_traced_end=float(tip_dict.get('accumulated_start_traced_len', 0.0) + tip_dict.get('accumulated_distances', 0.0)),
                                angle=float(tip_dict.get('angle', 0.0)),
                            )
                        ))
                        tip_dict['accumulated_points'].clear()
                        tip_dict['accumulated_distances'] = 0.0
                    
                    if depth < self.dendrite.depth-1:
                        self._spawn_children(tip_dict, new_tips)
                self.active_tips.remove(tip_dict); continue

            # Calculate distance and thickness values needed for collision detection
            dist_to_start_of_depth_level = sum(self.dendrite.branch_lengths[:depth]) if depth > 0 else 0
            safe_D_draw = self.dendrite.D if self.dendrite.D > 0 else 1.0
            current_dist_from_soma = dist_to_start_of_depth_level + tip_dict['traced_len']
            current_prop = current_dist_from_soma / self.dendrite.total_length
            current_thickness = max(1, self.dendrite.initial_thick * np.clip(1 - current_prop, 0, 1)**(1/safe_D_draw))

            # Determine growth direction using trajectory planning or traditional randomness,
            # with bounded retries for self-avoidance.
            p0_current_tip = tip_dict['pt']
            ang = tip_dict['angle']
            p1_potential_end = None
            max_angle_tries = 1 + (self.self_avoidance_retries if self.self_avoidance else 0)
            accepted_candidate = False

            for attempt_idx in range(max_angle_tries):
                self.last_step_stats['self_avoidance_attempts'] += 1
                t_dir = time.perf_counter()
                if self.trajectory_planning:
                    target_point, trajectory_angle = self._get_next_trajectory_target(tip_dict)
                    candidate_ang = trajectory_angle
                    if attempt_idx == 0:
                        tip_dict['last_replan_distance'] += current_growth_step
                    if not self.dendrite.det:
                        if self.dendrite.weave_type == 'Gauss':
                            candidate_ang += np.random.normal(0, self.dendrite.randomness * self.dendrite.curviness_mag * 0.05)
                        else:
                            candidate_ang += np.random.uniform(
                                -self.dendrite.randomness * self.dendrite.curviness_mag * 0.05,
                                self.dendrite.randomness * self.dendrite.curviness_mag * 0.05
                            )
                else:
                    candidate_ang = tip_dict['angle']
                    if not self.dendrite.det:
                        if self.dendrite.weave_type == 'Gauss':
                            candidate_ang += np.random.normal(0, self.dendrite.randomness * self.dendrite.curviness_mag * 0.1)
                        else:
                            candidate_ang += np.random.uniform(
                                -self.dendrite.randomness * self.dendrite.curviness_mag * 0.1,
                                self.dendrite.randomness * self.dendrite.curviness_mag * 0.1
                            )

                candidate_end = p0_current_tip + current_growth_step * np.array([np.cos(candidate_ang), np.sin(candidate_ang)])

                # Hard boundary reject for candidate direction.
                if (candidate_end[0] < 0 or candidate_end[0] >= self.network.width or
                    candidate_end[1] < 0 or candidate_end[1] >= self.network.height):
                    self.last_step_stats['growth_direction_select_s'] += time.perf_counter() - t_dir
                    continue
                self.last_step_stats['growth_direction_select_s'] += time.perf_counter() - t_dir

                if self.self_avoidance:
                    t_sa = time.perf_counter()
                    growth_vec = candidate_end - p0_current_tip
                    rr_try, cc_try = draw_line_pixels(
                        int(round(p0_current_tip[1])), int(round(p0_current_tip[0])),
                        int(round(candidate_end[1])), int(round(candidate_end[0]))
                    )
                    self_conflict = False
                    for r_try, c_try in zip(rr_try, cc_try):
                        if self._self_avoidance_conflict(
                            (c_try, r_try),
                            tip_dict,
                            current_thickness,
                            current_dist_from_soma,
                            growth_vec,
                        ):
                            self_conflict = True
                            break
                    if self_conflict:
                        self.last_step_stats['growth_self_avoidance_s'] += time.perf_counter() - t_sa
                        self.last_step_stats['self_avoidance_rejections'] += 1
                        continue
                    self.last_step_stats['growth_self_avoidance_s'] += time.perf_counter() - t_sa

                ang = candidate_ang
                p1_potential_end = candidate_end
                accepted_candidate = True
                break

            if not accepted_candidate:
                # Could not find a self-avoiding candidate direction; stop this tip.
                self.active_tips.remove(tip_dict)
                continue
            
            # Pixelized centerline for collision scan
            rr,cc = draw_line_pixels(int(round(p0_current_tip[1])),int(round(p0_current_tip[0])), \
                                     int(round(p1_potential_end[1])),int(round(p1_potential_end[0])))
            
            actual_end_point_this_step = p1_potential_end
            stopped_by_collision=False
            
            # Check collision along the centerline pixels
            # Use persistent sets stored in the tip to remember pass-through decisions
            encountered_dendrites = tip_dict['encountered_dendrites']
            pass_through_dendrites = tip_dict['pass_through_dendrites']
            if len(rr) > 0:
                t_coll = time.perf_counter()
                for r_px,c_px in zip(rr,cc):
                    # First check what's at this pixel without pass-through logic
                    collision_info = self._check_pixel_collision((c_px,r_px), current_thickness)
                    
                    if collision_info is not None:
                        collision_type, colliding_id, skeleton_contact = collision_info
                        idx_collision = list(zip(rr,cc)).index((r_px,c_px))
                        
                        if collision_type == 'boundary':
                            # Boundary collision always stops immediately
                            if idx_collision > 0:
                                actual_end_point_this_step = np.array([cc[idx_collision-1], rr[idx_collision-1]], float)
                            else:
                                actual_end_point_this_step = p0_current_tip
                            stopped_by_collision=True
                            break
                            
                        elif collision_type == 'soma':
                            # Soma collision - only stop when we actually touch the soma (skeleton_contact=True)
                            if skeleton_contact:
                                # Record soma-dendrite connection point
                                actual_end_point_this_step = np.array([c_px, r_px], float)
                                self.soma_connections.append((float(c_px), float(r_px), colliding_id, current_step))
                                stopped_by_collision = True
                                break
                            # If skeleton_contact=False, we're just passing through soma filled area, continue growing
                            
                        elif collision_type == 'dendrite':
                            # Dendrite collision - only process interactions when skeletons touch
                            if skeleton_contact:
                                # Check if we've already decided about this dendrite
                                if colliding_id not in encountered_dendrites:
                                    encountered_dendrites.add(colliding_id)
                                    
                                    # Three-stage decision process:
                                    # Stage 1: Can both dendrites coexist? (pass-through probability)
                                    # Stage 2: Should they fuse? (fusion probability + distance similarity)
                                    # Stage 3: If not, which one wins? (distance-based hierarchy)
                                    
                                    should_pass_through = False
                                    should_fuse = False
                                    
                                    # Stage 1: Roll for pass-through (can both dendrites coexist?)
                                    if not self.dendrite.det and random.random() < self.pass_prob:
                                        # Pass-through succeeds - both dendrites can coexist
                                        should_pass_through = True
                                    else:
                                        # Pass-through fails - check for fusion or competition
                                        current_dist_from_soma = dist_to_start_of_depth_level + tip_dict['traced_len']
                                        
                                        # Find the other neuron and estimate its distance at collision point
                                        other_neuron = None
                                        for n in self.network.neurons:
                                            if n.id_int == colliding_id:
                                                other_neuron = n
                                                break
                                        
                                        if other_neuron:
                                            # Calculate approximate distance from other neuron's soma to collision point
                                            other_soma_x, other_soma_y = other_neuron.position[:2]
                                            collision_dist_from_other_soma = np.sqrt((c_px - other_soma_x)**2 + (r_px - other_soma_y)**2)
                                            
                                            # Stage 2: Check for fusion (similar distances from respective somas)
                                            distance_diff = abs(current_dist_from_soma - collision_dist_from_other_soma)
                                            if (distance_diff <= self.fusion_distance_threshold and 
                                                not self.dendrite.det and random.random() < self.fusion_prob):
                                                # Fusion succeeds - both dendrites merge
                                                should_fuse = True
                                            else:
                                                # Stage 3: Use distance-based hierarchy for competition
                                                # Current dendrite wins if it's closer to its own soma
                                                if current_dist_from_soma < collision_dist_from_other_soma:
                                                    should_pass_through = True
                                                # If distances are very close, fall back to neuron ID hierarchy
                                                elif abs(current_dist_from_soma - collision_dist_from_other_soma) < 5.0:
                                                    should_pass_through = self.id_int < colliding_id
                                                # If current dendrite is farther, it loses (should_pass_through remains False)
                                    
                                    if should_fuse:
                                        # Record fusion point and stop this dendrite (it becomes part of fused dendrite)
                                        actual_end_point_this_step = np.array([c_px, r_px], float)
                                        self.fusion_points.append((float(c_px), float(r_px), colliding_id, current_step))
                                        stopped_by_collision = True
                                        break
                                    elif should_pass_through:
                                        # Either pass-through succeeded OR this dendrite won the competition
                                        pass_through_dendrites.add(colliding_id)
                                        continue
                                    else:
                                        # This dendrite lost the competition - stop and connect
                                        actual_end_point_this_step = np.array([c_px, r_px], float)
                                        self.connections.append((float(c_px), float(r_px), current_step))
                                        stopped_by_collision=True
                                        break
                                elif colliding_id in pass_through_dendrites:
                                    # We've already decided to pass through this dendrite
                                    continue
                                else:
                                    # We've encountered this dendrite before but it's not in pass_through list
                                    # This means we decided to stop at it before
                                    actual_end_point_this_step = np.array([c_px, r_px], float)
                                    self.connections.append((float(c_px), float(r_px), current_step))
                                    stopped_by_collision=True
                                    break
                            # If skeleton_contact=False, we're just passing through filled area, continue growing
                self.last_step_stats['growth_collision_scan_s'] += time.perf_counter() - t_coll
            
            # Calculate actual distance grown this step
            actual_grown_dist_this_step = np.linalg.norm(actual_end_point_this_step - p0_current_tip)

            # Calculate thickness for drawing (reuse already calculated values)
            dist_from_soma_to_p0 = current_dist_from_soma
            prop_s_draw = dist_from_soma_to_p0 / self.dendrite.total_length
            prop_e_draw = (dist_from_soma_to_p0 + actual_grown_dist_this_step) / self.dendrite.total_length
            
            thick_s_draw = max(1, self.dendrite.initial_thick * np.clip(1 - prop_s_draw,0,1)**(1/safe_D_draw))
            thick_e_draw = max(1, self.dendrite.initial_thick * np.clip(1 - prop_e_draw,0,1)**(1/safe_D_draw))
            avg_thick_draw = (thick_s_draw + thick_e_draw) / 2.0

            # Accumulate points for spline generation
            if actual_grown_dist_this_step > 1e-4:
                # Add start point if this is the first point in accumulation
                if len(tip_dict['accumulated_points']) == 0:
                    tip_dict['accumulated_points'].append(p0_current_tip.copy())
                    tip_dict['accumulated_start_traced_len'] = tip_dict['traced_len']
                
                # Add end point
                tip_dict['accumulated_points'].append(actual_end_point_this_step.copy())
                tip_dict['accumulated_distances'] += actual_grown_dist_this_step
                tip_dict['avg_accumulated_thickness'] = avg_thick_draw
                
                # Check if we should draw the accumulated splined segment
                should_draw_spline = (
                    tip_dict['accumulated_distances'] >= self.spline_segment_length or
                    stopped_by_collision or
                    tip_dict['traced_len'] + actual_grown_dist_this_step >= seg_len_for_this_depth - 1e-3
                )
                
                if should_draw_spline and len(tip_dict['accumulated_points']) >= 2:
                    chunk_start_abs = dist_to_start_of_depth_level + tip_dict.get('accumulated_start_traced_len', 0.0)
                    chunk_end_abs = chunk_start_abs + tip_dict['accumulated_distances']
                    splined_segments_to_draw.append((
                        tip_dict['accumulated_points'].copy(),
                        avg_thick_draw,
                        float(chunk_start_abs),
                        float(chunk_end_abs),
                        dict(
                            tip_id=tip_dict['id'],
                            parent_tip_id=tip_dict.get('parent_tip_id'),
                            parent_branch_path_dist=float(tip_dict.get('parent_branch_path_dist', 0.0)),
                            depth=depth,
                            tip_traced_start=float(tip_dict.get('accumulated_start_traced_len', 0.0)),
                            tip_traced_end=float(tip_dict.get('accumulated_start_traced_len', 0.0) + tip_dict.get('accumulated_distances', 0.0)),
                            angle=float(ang),
                        )
                    ))
                    # Reset accumulation but keep last point as start of next segment
                    last_point = tip_dict['accumulated_points'][-1].copy()
                    tip_dict['accumulated_points'] = [last_point]
                    tip_dict['accumulated_distances'] = 0.0
                    tip_dict['accumulated_start_traced_len'] = tip_dict['traced_len']
            
            tip_dict['pt'] = actual_end_point_this_step
            tip_dict['traced_len'] += actual_grown_dist_this_step
            
            # Update angle for next iteration - use the calculated angle from trajectory planning
            tip_dict['angle'] = ang
            
            # Check for segment completion or collision
            if stopped_by_collision or tip_dict['traced_len'] >= seg_len_for_this_depth - 1e-3:
                # Draw any remaining accumulated points
                if len(tip_dict['accumulated_points']) > 1:
                    chunk_start_abs = dist_to_start_of_depth_level + tip_dict.get('accumulated_start_traced_len', 0.0)
                    chunk_end_abs = chunk_start_abs + tip_dict['accumulated_distances']
                    splined_segments_to_draw.append((
                        tip_dict['accumulated_points'].copy(),
                        avg_thick_draw,
                        float(chunk_start_abs),
                        float(chunk_end_abs),
                        dict(
                            tip_id=tip_dict['id'],
                            parent_tip_id=tip_dict.get('parent_tip_id'),
                            parent_branch_path_dist=float(tip_dict.get('parent_branch_path_dist', 0.0)),
                            depth=depth,
                            tip_traced_start=float(tip_dict.get('accumulated_start_traced_len', 0.0)),
                            tip_traced_end=float(tip_dict.get('accumulated_start_traced_len', 0.0) + tip_dict.get('accumulated_distances', 0.0)),
                            angle=float(ang),
                        )
                    ))
                
                self.active_tips.remove(tip_dict)
                
                if depth < self.dendrite.depth-1 and not stopped_by_collision:
                    self._spawn_children(tip_dict, new_tips, parent_final_angle=ang)
        
        # Draw all splined segments after collision detection is complete
        t_draw = time.perf_counter()
        for accumulated_points, thickness, path_start_dist, path_end_dist, seg_meta in splined_segments_to_draw:
            self._draw_splined_segment(
                accumulated_points,
                thickness,
                path_start_dist=path_start_dist,
                path_end_dist=path_end_dist,
                seg_meta=seg_meta
            )
        self.last_step_stats['growth_draw_record_s'] += time.perf_counter() - t_draw
            
        self.active_tips.extend(new_tips)
        if not self.active_tips:
            self.is_growing=False
            
    def _spawn_children(self, parent_tip_dict, new_tips_list, parent_final_angle=None):
        """ Helper to spawn child tips """
        depth = parent_tip_dict['depth']
        p1 = parent_tip_dict['pt']
        
        base_angle_for_children = parent_final_angle if parent_final_angle is not None else parent_tip_dict['angle']

        k = self.dendrite.mean_branches
        if not self.dendrite.det and self.dendrite.mean_branches > 0:
            k = max(1,int(np.round(np.random.normal(self.dendrite.mean_branches,1))))
        
        for i in range(k):
            new_angle = base_angle_for_children + self.dendrite.branch_angle*(i-(k-1)/2.0)
            if not self.dendrite.det:
                 new_angle += np.random.normal(0, self.dendrite.randomness * 0.2)

            new_tips_list.append(dict(id=self._tip_id, pt=p1.copy(), angle=new_angle, depth=depth+1, traced_len=0.0,
                                      parent_tip_id=parent_tip_dict.get('id'),
                                      parent_branch_path_dist=float(sum(self.dendrite.branch_lengths[:depth]) + parent_tip_dict.get('traced_len', 0.0)),
                                      encountered_dendrites=set(), pass_through_dendrites=set(),
                                      accumulated_points=[], accumulated_distances=0.0,
                                      accumulated_start_traced_len=0.0,
                                      # Trajectory planning data
                                      trajectory_points=None,  # List of target points along planned trajectory
                                      trajectory_index=0,      # Current position in trajectory
                                      last_replan_distance=0.0  # Distance since last trajectory replan
                                      ))
            self._tip_id+=1

    def _draw_splined_segment(self, accumulated_points: List[np.ndarray], avg_thickness: float,
                              path_start_dist: float=None, path_end_dist: float=None,
                              seg_meta: Dict=None):
        """
        Take accumulated points and draw a splined curve through them.
        """
        if len(accumulated_points) < 2:
            return
            
        if not self.use_splines or self.dendrite.det:
            # Draw straight segments if splines disabled or deterministic
            self._record_branch_segment(accumulated_points, avg_thickness, path_start_dist, path_end_dist, seg_meta=seg_meta)
            for i in range(len(accumulated_points) - 1):
                self._draw_segment(accumulated_points[i], accumulated_points[i+1], avg_thickness)
            return
            
        # Calculate total length of accumulated path
        total_length = 0
        for i in range(len(accumulated_points) - 1):
            total_length += np.linalg.norm(accumulated_points[i+1] - accumulated_points[i])
            
        if total_length < 2.0:  # Too short for meaningful spline
            self._record_branch_segment(accumulated_points, avg_thickness, path_start_dist, path_end_dist, seg_meta=seg_meta)
            for i in range(len(accumulated_points) - 1):
                self._draw_segment(accumulated_points[i], accumulated_points[i+1], avg_thickness)
            return
            
        # Use spline to create smooth curve
        start_pt = accumulated_points[0]
        end_pt = accumulated_points[-1]
        
        try:
            # Generate splined curve
            xs, ys = self.dendrite.intra_branch_weave(
                start_pt[0], start_pt[1], end_pt[0], end_pt[1], total_length
            )
            drawn_points = [np.array([x, y], dtype=float) for x, y in zip(xs, ys)]
            self._record_branch_segment(drawn_points, avg_thickness, path_start_dist, path_end_dist, seg_meta=seg_meta)
            
            # Draw the splined curve as segments
            for i in range(len(xs) - 1):
                p0 = np.array([xs[i], ys[i]])
                p1 = np.array([xs[i+1], ys[i+1]])
                self._draw_segment(p0, p1, avg_thickness)
                
        except Exception:
            self._record_branch_segment(accumulated_points, avg_thickness, path_start_dist, path_end_dist, seg_meta=seg_meta)
            # Fallback to straight segments if spline fails
            for i in range(len(accumulated_points) - 1):
                self._draw_segment(accumulated_points[i], accumulated_points[i+1], avg_thickness)

    def _record_branch_segment(self, points: List[np.ndarray], thickness: float,
                               path_start_dist: float=None, path_end_dist: float=None,
                               seg_meta: Dict=None):
        """Store drawable segment history so masks can be rebuilt/pruned in place."""
        if points is None or len(points) < 2:
            return
        pts = [np.array(p, dtype=float).copy() for p in points]

        # Estimate path distance span if not provided.
        if path_start_dist is None or path_end_dist is None:
            seg_len = 0.0
            for i in range(len(pts) - 1):
                seg_len += float(np.linalg.norm(pts[i+1] - pts[i]))
            path_start_dist = float(self.max_recorded_path_length)
            path_end_dist = float(path_start_dist + seg_len)

        path_start_dist = float(path_start_dist)
        path_end_dist = float(max(path_start_dist, path_end_dist))
        self.max_recorded_path_length = max(self.max_recorded_path_length, path_end_dist)

        seg_record = {
            'points': pts,
            'thickness': float(max(1.0, thickness)),
            'path_start_dist': path_start_dist,
            'path_end_dist': path_end_dist,
            'segment_index': len(self.branch_data),
        }
        if seg_meta:
            seg_record.update(dict(seg_meta))
        self.branch_data.append(seg_record)
        if getattr(self.network, 'incremental_branch_mode', False):
            # In the online simulation path, newly recorded segments are immediately visible.
            self.network.incremental_add_segment_record(self.id_int, seg_record)

    def _clip_polyline_fraction_local(self, points: List[np.ndarray], fraction: float) -> List[np.ndarray]:
        """Local polyline clip helper used for prune-and-forget frontier extraction."""
        if not points or len(points) < 2:
            return []
        fraction = float(np.clip(fraction, 0.0, 1.0))
        if fraction <= 0.0:
            return []
        if fraction >= 1.0:
            return [np.array(p, dtype=float).copy() for p in points]

        pts = [np.array(p, dtype=float) for p in points]
        seg_lengths = [float(np.linalg.norm(pts[i+1] - pts[i])) for i in range(len(pts) - 1)]
        total_len = float(sum(seg_lengths))
        if total_len <= 1e-9:
            return [pts[0].copy(), pts[-1].copy()]

        target = total_len * fraction
        acc = 0.0
        clipped = [pts[0].copy()]
        for i, seg_len in enumerate(seg_lengths):
            p0 = pts[i]
            p1 = pts[i + 1]
            if acc + seg_len < target - 1e-9:
                clipped.append(p1.copy())
                acc += seg_len
                continue
            rem = max(0.0, target - acc)
            if seg_len <= 1e-9:
                clipped.append(p1.copy())
            else:
                alpha = float(np.clip(rem / seg_len, 0.0, 1.0))
                clipped.append((1.0 - alpha) * p0 + alpha * p1)
            break
        if len(clipped) < 2:
            clipped.append(pts[1].copy())
        return clipped

    def _make_frontier_tip(self, seg: Dict) -> Dict:
        points = seg.get('points', [])
        if not points or len(points) < 2:
            return None
        p_last = np.array(points[-1], dtype=float)
        p_prev = np.array(points[-2], dtype=float)
        d = p_last - p_prev
        if np.linalg.norm(d) <= 1e-9:
            angle = float(seg.get('angle', 0.0))
        else:
            angle = float(np.arctan2(d[1], d[0]))

        return dict(
            id=int(seg.get('tip_id', self._tip_id)),
            pt=p_last.copy(),
            angle=angle,
            depth=int(seg.get('depth', 0)),
            traced_len=float(seg.get('tip_traced_end', 0.0)),
            parent_tip_id=seg.get('parent_tip_id'),
            parent_branch_path_dist=float(seg.get('parent_branch_path_dist', 0.0)),
            encountered_dendrites=set(),
            pass_through_dendrites=set(),
            accumulated_points=[],
            accumulated_distances=0.0,
            accumulated_start_traced_len=0.0,
            trajectory_points=None,
            trajectory_index=0,
            last_replan_distance=0.0,
        )

    def prune_and_forget(self, cutoff_length: float, eps: float=1e-6):
        """
        Permanently delete recorded branch history beyond cutoff_length and reactivate
        growth tips at the prune frontier so future regrowth can follow new paths.
        """
        cutoff = float(max(0.0, cutoff_length))

        current_max = float(getattr(self, 'max_recorded_path_length', 0.0))
        if cutoff >= current_max - eps:
            # No truncation needed.
            return

        self.set_visible_total_length(cutoff)

        if not self.branch_data:
            # No recorded geometry yet: keep current tips and re-sprout any lost primaries.
            self.restore_missing_primary_tips()
            return

        prev_active_tips = list(self.active_tips)
        incremental_mode = bool(getattr(self.network, 'incremental_branch_mode', False))
        old_history = self.branch_data

        # Build a tip lineage map and prune tips by ancestry/branchpoint, then prune each
        # tip as a contiguous prefix. This avoids leaving orphan chunks behind.
        tip_segments = {}
        tip_meta = {}
        for seg in self.branch_data:
            tip_id = seg.get('tip_id', None)
            if tip_id is None:
                # Legacy/no-lineage segment: treat as its own root tip bucket.
                tip_id = f"legacy_{seg.get('segment_index', len(tip_segments))}"
            tip_segments.setdefault(tip_id, []).append(seg)
            if tip_id not in tip_meta:
                tip_meta[tip_id] = dict(
                    parent_tip_id=seg.get('parent_tip_id', None),
                    parent_branch_path_dist=float(seg.get('parent_branch_path_dist', 0.0)),
                    depth=int(seg.get('depth', 0)),
                )

        # Topologically evaluate which tips are still reachable after prune.
        tip_keepable = {}
        unresolved = set(tip_segments.keys())
        while unresolved:
            progressed = False
            for tip_id in list(unresolved):
                meta = tip_meta.get(tip_id, {})
                parent_tip_id = meta.get('parent_tip_id', None)
                branch_dist = float(meta.get('parent_branch_path_dist', 0.0))
                if parent_tip_id is None or parent_tip_id == tip_id:
                    tip_keepable[tip_id] = True
                    unresolved.remove(tip_id)
                    progressed = True
                    continue
                if parent_tip_id not in tip_segments:
                    tip_keepable[tip_id] = cutoff >= branch_dist - eps
                    unresolved.remove(tip_id)
                    progressed = True
                    continue
                if parent_tip_id in tip_keepable:
                    tip_keepable[tip_id] = bool(tip_keepable[parent_tip_id]) and (cutoff >= branch_dist - eps)
                    unresolved.remove(tip_id)
                    progressed = True
            if not progressed:
                # Cycle/corruption fallback: mark remaining tips non-keepable.
                for tip_id in unresolved:
                    tip_keepable[tip_id] = False
                break

        new_history = []
        last_by_tip = {}
        truncated_tips = set()   # tips that lost at least part of their recorded path

        for tip_id, segs in tip_segments.items():
            if not tip_keepable.get(tip_id, False):
                continue

            segs_sorted = sorted(
                segs,
                key=lambda s: (
                    float(s.get('tip_traced_start', 0.0)),
                    float(s.get('path_start_dist', 0.0)),
                    int(s.get('segment_index', 0)),
                )
            )

            # Keep a contiguous prefix only.
            for seg in segs_sorted:
                path_start = float(seg.get('path_start_dist', 0.0))
                path_end = float(seg.get('path_end_dist', path_start))
                if cutoff <= path_start + eps:
                    truncated_tips.add(tip_id)
                    break

                kept_seg = seg   # unchanged segment: reuse object + cached rasterisation
                clipped = False
                if cutoff < path_end - eps:
                    frac = (cutoff - path_start) / max(eps, (path_end - path_start))
                    clipped_points = self._clip_polyline_fraction_local(seg.get('points', []), frac)
                    if len(clipped_points) < 2:
                        truncated_tips.add(tip_id)
                        break
                    kept_seg = dict(seg)
                    kept_seg.pop('filled_idx', None)
                    kept_seg.pop('skeleton_idx', None)
                    kept_seg['points'] = clipped_points
                    kept_seg['path_end_dist'] = cutoff
                    clipped = True
                    truncated_tips.add(tip_id)

                    if 'tip_traced_start' in kept_seg and 'tip_traced_end' in kept_seg:
                        ts = float(kept_seg['tip_traced_start'])
                        te = float(kept_seg['tip_traced_end'])
                        kept_seg['tip_traced_end'] = ts + frac * (te - ts)

                new_history.append(kept_seg)
                last_by_tip[tip_id] = kept_seg

                if clipped:
                    break

        self.branch_data = new_history
        self.max_recorded_path_length = max(
            [float(seg.get('path_end_dist', 0.0)) for seg in self.branch_data] + [0.0]
        )
        if incremental_mode:
            kept_ids = {id(seg) for seg in self.branch_data}
            old_ids = {id(seg) for seg in old_history}
            removed = [seg for seg in old_history if id(seg) not in kept_ids]
            added = [seg for seg in self.branch_data if id(seg) not in old_ids]
            self.network.incremental_batch_update(self.id_int, removed, added)

        frontier_tips = []
        for seg in last_by_tip.values():
            seg_end = float(seg.get('path_end_dist', 0.0))
            if abs(seg_end - cutoff) > 1e-3:
                continue
            depth = int(seg.get('depth', 0))
            if depth >= self.dendrite.depth:
                continue
            tip = self._make_frontier_tip(seg)
            if tip is not None:
                frontier_tips.append(tip)

        # Keep tips that were still growing and are untouched by the prune (shorter
        # branches of the same neuron); previously these were silently dropped.
        dedup = {}
        for tip in prev_active_tips:
            tid = tip['id']
            if tid in tip_segments:
                keep = tip_keepable.get(tid, False) and tid not in truncated_tips
            else:
                # No recorded geometry yet: keep if its parent branch point survives.
                parent = tip.get('parent_tip_id')
                keep = (parent is None) or (
                    tip_keepable.get(parent, parent not in tip_segments)
                    and cutoff >= float(tip.get('parent_branch_path_dist', 0.0)) - eps
                    and (parent not in truncated_tips
                         or float(tip.get('parent_branch_path_dist', 0.0)) <= cutoff - eps))
            if keep:
                dedup[tid] = tip
        # Frontier tips at the prune cut (new growth cones on truncated branches).
        for tip in frontier_tips:
            dedup.setdefault(tip['id'], tip)
        self.active_tips = list(dedup.values())

        # Ensure future spawned tips continue using unique IDs.
        if self.active_tips:
            self._tip_id = max(self._tip_id, max(int(t['id']) for t in self.active_tips) + 1)
        self.restore_missing_primary_tips()

    def _new_root_tip(self, tid: int, pt: np.ndarray, angle: float) -> Dict:
        return dict(id=tid, pt=np.array(pt, dtype=float).copy(), angle=float(angle), depth=0, traced_len=0.0,
                    parent_tip_id=None, parent_branch_path_dist=0.0,
                    encountered_dendrites=set(), pass_through_dendrites=set(),
                    accumulated_points=[], accumulated_distances=0.0,
                    accumulated_start_traced_len=0.0,
                    trajectory_points=None, trajectory_index=0,
                    last_replan_distance=0.0)

    def sprout_new_primary(self, n_candidates: int = 24, probe_px: int = 20) -> bool:
        """
        Grow a brand-new primary dendrite from the soma, in the least-crowded direction.
        Used when the neuron is below its homeostatic target but every existing tip has been
        stopped (collisions / canvas edge), so it would otherwise be permanently under-driven.
        Returns True if a tip was added.
        """
        if self.dendrite.depth <= 0:
            return False
        xs, ys = np.asarray(self.dendrite.soma.x, float), np.asarray(self.dendrite.soma.y, float)
        cx, cy = self.position[:2]
        mask = self.network.filled_mask
        H, W = mask.shape
        contour_angs = (np.arctan2(ys - cy, xs - cx) + 2 * np.pi) % (2 * np.pi)
        a0 = np.random.uniform(0, 2 * np.pi)
        best = None
        best_occ = np.inf
        d = np.arange(3.0, 3.0 + probe_px, 2.0)
        for k in range(n_candidates):
            a = (a0 + 2 * np.pi * k / n_candidates) % (2 * np.pi)
            idx = int(np.argmin(np.abs((contour_angs - a + np.pi) % (2 * np.pi) - np.pi)))
            px = np.round(xs[idx] + d * np.cos(a)).astype(int)
            py = np.round(ys[idx] + d * np.sin(a)).astype(int)
            inb = (px >= 0) & (px < W) & (py >= 0) & (py < H)
            if inb.sum() < len(d):
                continue  # would leave the canvas almost immediately
            occ = int(np.count_nonzero(mask[py, px]))
            if occ < best_occ:
                best_occ, best = occ, (np.array([xs[idx], ys[idx]], float), a)
        if best is None:
            return False
        pt, a = best
        if not hasattr(self, '_primary_tip_templates'):
            self._primary_tip_templates = {}
        tid = self._tip_id
        self._tip_id += 1
        self._primary_tip_templates[tid] = (pt.copy(), float(a))
        if not hasattr(self, '_resprout_ids'):
            self._resprout_ids = set()
        self._resprout_ids.add(tid)   # re-sprouted primaries are NOT restored after pruning
        self.active_tips.append(self._new_root_tip(tid, pt, a))
        self.n_resprouts = getattr(self, 'n_resprouts', 0) + 1
        self.is_growing = True
        return True

    def retract_blocked_resprouts(self, min_length: float) -> int:
        """
        Remove re-sprouted primaries that stopped growing (all tips blocked) while still
        shorter than `min_length` px, so failed sprouts don't remain as stubs on the soma.
        Returns the number of sprouts retracted.
        """
        resprouts = getattr(self, '_resprout_ids', None)
        if not resprouts:
            return 0
        active_ids = {t['id'] for t in self.active_tips}
        active_parents = {t.get('parent_tip_id') for t in self.active_tips}
        retracted = 0
        for tid in list(resprouts):
            if tid in active_ids or tid in active_parents:
                continue                      # still growing
            segs = [seg for seg in self.branch_data if seg.get('tip_id') == tid]
            if any(seg.get('parent_tip_id') == tid for seg in self.branch_data):
                continue                      # it branched -> a real dendrite, keep it
            length = max([float(seg.get('path_end_dist', 0.0)) for seg in segs] + [0.0])
            if length >= min_length:
                continue
            if segs:
                if getattr(self.network, 'incremental_branch_mode', False):
                    self.network.incremental_batch_update(self.id_int, segs, [])
                drop = {id(seg) for seg in segs}
                self.branch_data = [seg for seg in self.branch_data if id(seg) not in drop]
            resprouts.discard(tid)
            self._primary_tip_templates.pop(tid, None)
            retracted += 1
        if retracted:
            self.max_recorded_path_length = max(
                [float(seg.get('path_end_dist', 0.0)) for seg in self.branch_data] + [0.0])
        return retracted

    def restore_missing_primary_tips(self, angle_jitter: float = np.pi / 12) -> int:
        """
        Re-sprout primary dendrites that have been pruned all the way back to the soma
        (no recorded geometry and no active tip left). Without this, a neuron pruned to a
        bare soma could never regrow. Returns the number of tips restored.
        """
        templates = getattr(self, '_primary_tip_templates', {})
        if not templates or self.dendrite.depth <= 0:
            self.is_growing = bool(self.active_tips)
            return 0
        recorded_roots = {seg.get('tip_id') for seg in self.branch_data}
        active_ids = {t['id'] for t in self.active_tips}
        resprouts = getattr(self, '_resprout_ids', set())
        restored = 0
        for tid, (pt, a) in list(templates.items()):
            if tid in recorded_roots or tid in active_ids:
                continue
            if tid in resprouts:
                # A re-sprouted primary that has been pruned away entirely is gone for good
                # (restoring it would leave the soma covered in stubs).
                resprouts.discard(tid)
                del templates[tid]
                continue
            if not self.dendrite.det:
                a = a + np.random.uniform(-angle_jitter, angle_jitter)
            self.active_tips.append(self._new_root_tip(tid, pt, a))
            restored += 1
        self.is_growing = bool(self.active_tips)
        return restored

################################################################################
#                               NETWORK CLASS                                  #
################################################################################

class Network:
    def __init__(self, width:int, height:int, num_neurons:int, neuron_params:dict, edge_margin:int=50,
                 sim_props:dict=None, deterministic=False, seed_coordinates=None, network_id=None, simage_props=None,
                 network_params=None):
        self.width=width; self.height=height; self.edge_margin=edge_margin
        self.neuron_params=neuron_params
        self.deterministic=deterministic
        self.neurons: List[Neuron]=[]
        self.network_id = network_id if network_id is not None else "network"
        
        # Extract placement parameters from network_params
        network_params = network_params or {}
        self.manual_buffer_distance = network_params.get('manual_buffer_distance', None)  # Override computed buffer distance
        self.jitter_factor = network_params.get('jitter_factor', 0.2)  # How much to randomize positions (0=no jitter, 1=max jitter)

        # raster layers
        self.filled_mask   = np.zeros((height,width), np.uint16)
        self.skeleton_mask = np.zeros((height,width), np.uint16)
        self.somas_mask    = np.zeros((height,width), np.uint16)
        self.soma_skeleton_mask = np.zeros((height,width), np.uint16)
        self.simage        = np.zeros((height,width), np.uint8)
        self.outline_mask  = np.zeros((height,width), np.uint8)
        self.last_rebuild_history_stats = {
            'segments_replayed': 0,
            'segments_clipped': 0,
            'segments_skipped': 0,
        }
        self.incremental_branch_mode = True
        self._branch_filled_counts = [np.zeros((height, width), np.int32) for _ in range(num_neurons)]
        self._branch_skeleton_counts = [np.zeros((height, width), np.int32) for _ in range(num_neurons)]

        # time controls
        sim_props = sim_props or {}
        self.dt   = sim_props.get('dt',1.0)
        self.max_steps = sim_props.get('max_growth_steps',1000)
        self.max_growth_steps_safety_cap = sim_props.get('max_growth_steps_safety_cap', None)
        self.record_history_every_n_steps = sim_props.get('record_history_every_n_steps',0)
        self.growth_history_masks = []
        
        # Morphological smoothing parameters
        self.morphological_smoothing = sim_props.get('morphological_smoothing', False)
        self.local_smoothing_method = sim_props.get('local_smoothing_method', 'closing')  # 'closing', 'gaussian', 'opening', 'distance'
        self.smoothing_radius = sim_props.get('smoothing_radius', 3)
        self.smoothing_zone_radius = sim_props.get('smoothing_zone_radius', 10)  # Local area around connection/fusion points
        self.gaussian_sigma = sim_props.get('gaussian_sigma', 1.0)  # For Gaussian smoothing

        current_simage_props = simage_props if simage_props is not None else {}
        self.simage_blur_kernel_size = current_simage_props.get('simage_blur_kernel_size', (5,5))
        self.simage_blur_sigma = current_simage_props.get('simage_blur_sigma', 0)

        # Seeding
        self.seed_coordinates = seed_coordinates
        self._seed_neurons(num_neurons)

    # ------------------------------------------------------------- seeding
    def _seed_neurons(self, num_neurons:int):
        if self.seed_coordinates is not None:
            assert len(self.seed_coordinates) == num_neurons, \
                f"Number of manual coordinates ({len(self.seed_coordinates)}) must match num_neurons"
            for idx, seed_pos_tuple in enumerate(self.seed_coordinates):
                n = Neuron(seed_pos_tuple, self.neuron_params, self, f'{self.network_id}_N{idx+1}')
                n.id_int = idx + 1
                self.neurons.append(n)
                # Draw soma immediately
                f,s = n.dendrite.soma.binary_masks((self.height,self.width), n.id_int)
                self.somas_mask[f>0]=f[f>0]
                self.filled_mask[f>0]=f[f>0]
                self.skeleton_mask[s>0]=s[s>0]
                self.soma_skeleton_mask[s>0]=s[s>0]

        else:
            mean_r = self.neuron_params['mean_soma_radius']
            std_r  = self.neuron_params['std_soma_radius']
            
            # Use manual buffer distance if provided, otherwise compute it
            if self.manual_buffer_distance is not None:
                buf = self.manual_buffer_distance
                print(f"Using manual buffer distance: {buf}")
            else:
                buf = compute_max_buffer_distance(self.width,self.height,num_neurons,mean_r,std_r,self.edge_margin)
                print(f"Computed buffer distance: {buf}")
                
            cell_diameter = 2*(mean_r+2*std_r)+buf
            
            if cell_diameter <=0:
                raise ValueError("Calculated cell diameter for seeding is not positive. Check soma radii and buffer.")

            cols   = int((self.width -2*self.edge_margin)//cell_diameter)
            rows   = int((self.height-2*self.edge_margin)//cell_diameter)

            if cols == 0 or rows == 0 or num_neurons > cols*rows:
                raise ValueError(f"Cannot fit {num_neurons} in a {cols}x{rows} grid with cell size {cell_diameter:.2f}. Adjust params.")

            possible_positions = random.sample([(i,j) for i in range(cols) for j in range(rows)], num_neurons)
            for idx,(i,j) in enumerate(possible_positions):
                cell_center_x = self.edge_margin + (i+0.5)*cell_diameter
                cell_center_y = self.edge_margin + (j+0.5)*cell_diameter
                
                # Configurable jitter amount based on jitter_factor
                jitter_allowance_x = cell_diameter * self.jitter_factor
                jitter_allowance_y = cell_diameter * self.jitter_factor
                
                max_s_r = mean_r + 2*std_r
                
                pos_x = cell_center_x + random.uniform(-jitter_allowance_x/2, jitter_allowance_x/2)
                pos_y = cell_center_y + random.uniform(-jitter_allowance_y/2, jitter_allowance_y/2)

                pos_x = np.clip(pos_x, self.edge_margin + max_s_r, self.width - self.edge_margin - max_s_r)
                pos_y = np.clip(pos_y, self.edge_margin + max_s_r, self.height - self.edge_margin - max_s_r)

                n  = Neuron((pos_x,pos_y), self.neuron_params, self, f'{self.network_id}_N{idx+1}')
                n.id_int = idx+1
                self.neurons.append(n)
                # draw soma immediately
                f,s = n.dendrite.soma.binary_masks((self.height,self.width), n.id_int)
                self.somas_mask[f>0]=f[f>0]
                self.filled_mask[f>0]=f[f>0]
                self.skeleton_mask[s>0]=s[s>0]
                self.soma_skeleton_mask[s>0]=s[s>0]

    # ------------------------------------------------------------- growth loop
    def grow_network(self):
        def _record_growth_history(step_num:int):
            if not (self.record_history_every_n_steps > 0 and (step_num % self.record_history_every_n_steps == 0)):
                return
            # Collect connections and fusions that existed at this step
            current_connections = []
            current_fusions = []
            current_soma_connections = []
            for n in self.neurons:
                for conn in n.connections:
                    if conn[2] <= step_num:  # conn = (x, y, step)
                        current_connections.append(conn)
                for fusion in n.fusion_points:
                    if fusion[3] <= step_num:  # fusion = (x, y, other_id, step)
                        current_fusions.append(fusion)
                for soma_conn in n.soma_connections:
                    if soma_conn[3] <= step_num:  # soma_conn = (x, y, soma_id, step)
                        current_soma_connections.append(soma_conn)

            self.growth_history_masks.append({
                'step': step_num,
                'filled_mask': self.filled_mask.copy(),
                'skeleton_mask': self.skeleton_mask.copy(),
                'somas_mask': self.somas_mask.copy(),
                'connections': current_connections,
                'fusions': current_fusions,
                'soma_connections': current_soma_connections
            })

        active_neurons_this_step = False
        step_num = 0

        if self.max_steps is None:
            # Unbounded growth mode: stop only when all neurons stop growing.
            pbar = tqdm(total=None, desc="Growth step")
            while True:
                active_neurons_this_step = False
                for n in self.neurons:
                    if n.is_growing:
                        n.step(self.dt, current_step=step_num)
                        if n.is_growing:
                            active_neurons_this_step = True

                _record_growth_history(step_num)
                pbar.update(1)

                if not active_neurons_this_step:
                    print(f"Network growth stabilized at step {step_num}")
                    break

                if (self.max_growth_steps_safety_cap is not None and
                    step_num + 1 >= int(self.max_growth_steps_safety_cap)):
                    print(
                        f"Network growth hit safety cap ({self.max_growth_steps_safety_cap}) "
                        "before stabilization."
                    )
                    break

                step_num += 1
            pbar.close()
        else:
            for step_num in tqdm(range(self.max_steps), desc="Growth step"):
                active_neurons_this_step=False
                for n in self.neurons:
                    if n.is_growing:
                        n.step(self.dt, current_step=step_num)
                        if n.is_growing:
                            active_neurons_this_step=True

                _record_growth_history(step_num)

                if not active_neurons_this_step:
                    print(f"Network growth stabilized at step {step_num}")
                    break
            else:
                 if active_neurons_this_step:
                     print(f"Network growth reached max_growth_steps ({self.max_steps}) with some neurons still active.")
                 else:
                     print(f"Network growth stabilized at max_growth_steps ({self.max_steps}).")
        
        self._finalize_masks()

    def _finalize_masks(self):
        """Finalize additional masks like outline and simage after growth."""
        
        # Apply local morphological smoothing if enabled
        if self.morphological_smoothing:
            self._apply_local_smoothing()
        
        self.outline_mask = (sobel(self.filled_mask > 0) > 0).astype(np.uint8)
        if np.any(self.filled_mask):
            binary_filled_for_blur = (self.filled_mask > 0).astype(np.float32) * 255
            k_size_x = self.simage_blur_kernel_size[0]
            k_size_y = self.simage_blur_kernel_size[1]
            if k_size_x % 2 == 0: k_size_x +=1
            if k_size_y % 2 == 0: k_size_y +=1
            
            blurred_simage = cv2.GaussianBlur(binary_filled_for_blur, (k_size_x, k_size_y), self.simage_blur_sigma)
            if blurred_simage.max() > 0:
                blurred_simage = (blurred_simage / blurred_simage.max()) * 255
            self.simage = blurred_simage.astype(np.uint8)
        else:
            self.simage.fill(0)

    def get_recorded_branch_length_limits(self):
        """Return the maximum recorded path distance per neuron."""
        return np.array([float(getattr(n, 'max_recorded_path_length', 0.0)) for n in self.neurons], dtype=float)

    def _clip_polyline_to_fraction(self, points: List[np.ndarray], fraction: float) -> List[np.ndarray]:
        """Clip a polyline to a prefix fraction of its path length."""
        if not points or len(points) < 2:
            return []
        fraction = float(np.clip(fraction, 0.0, 1.0))
        if fraction <= 0.0:
            return []
        if fraction >= 1.0:
            return [np.array(p, dtype=float).copy() for p in points]

        pts = [np.array(p, dtype=float) for p in points]
        seg_lengths = [float(np.linalg.norm(pts[i+1] - pts[i])) for i in range(len(pts) - 1)]
        total_len = float(sum(seg_lengths))
        if total_len <= 1e-9:
            return [pts[0].copy(), pts[-1].copy()]

        target_len = fraction * total_len
        accumulated = 0.0
        clipped = [pts[0].copy()]

        for i, seg_len in enumerate(seg_lengths):
            p0 = pts[i]
            p1 = pts[i + 1]
            if accumulated + seg_len < target_len - 1e-9:
                clipped.append(p1.copy())
                accumulated += seg_len
                continue

            remaining = max(0.0, target_len - accumulated)
            if seg_len <= 1e-9:
                clipped.append(p1.copy())
            else:
                alpha = float(np.clip(remaining / seg_len, 0.0, 1.0))
                clipped.append((1.0 - alpha) * p0 + alpha * p1)
            break

        if len(clipped) < 2:
            clipped.append(pts[1].copy())
        return clipped

    def _draw_history_polyline(self, neuron_id: int, points: List[np.ndarray], thickness: float):
        """Draw a recorded polyline segment into filled/skeleton masks."""
        if points is None or len(points) < 2:
            return
        thickness_i = int(max(1, round(float(thickness))))
        for i in range(len(points) - 1):
            p0i = tuple(np.round(points[i]).astype(int))
            p1i = tuple(np.round(points[i+1]).astype(int))
            if (0 <= p0i[0] < self.width and 0 <= p0i[1] < self.height and
                0 <= p1i[0] < self.width and 0 <= p1i[1] < self.height):
                cv2.line(self.filled_mask, p0i, p1i, int(neuron_id), thickness=thickness_i)
                cv2.line(self.skeleton_mask, p0i, p1i, int(neuron_id), thickness=1)

    def _rasterize_segment_indices(self, points: List[np.ndarray], thickness: float):
        """Rasterize a polyline into compact global pixel index arrays (filled + skeleton)."""
        if points is None or len(points) < 2:
            return np.array([], dtype=np.int64), np.array([], dtype=np.int64)

        pts = np.array(points, dtype=float)
        th = int(max(1, round(float(thickness))))
        margin = th + 2
        x_min = max(0, int(np.floor(np.min(pts[:, 0]))) - margin)
        x_max = min(self.width - 1, int(np.ceil(np.max(pts[:, 0]))) + margin)
        y_min = max(0, int(np.floor(np.min(pts[:, 1]))) - margin)
        y_max = min(self.height - 1, int(np.ceil(np.max(pts[:, 1]))) + margin)
        if x_max < x_min or y_max < y_min:
            return np.array([], dtype=np.int64), np.array([], dtype=np.int64)

        h = y_max - y_min + 1
        w = x_max - x_min + 1
        local_filled = np.zeros((h, w), np.uint8)
        local_skel = np.zeros((h, w), np.uint8)

        for i in range(len(points) - 1):
            p0 = np.array(points[i], dtype=float)
            p1 = np.array(points[i + 1], dtype=float)
            p0i = (int(round(p0[0])) - x_min, int(round(p0[1])) - y_min)
            p1i = (int(round(p1[0])) - x_min, int(round(p1[1])) - y_min)
            if (0 <= p0i[0] < w and 0 <= p0i[1] < h and 0 <= p1i[0] < w and 0 <= p1i[1] < h):
                cv2.line(local_filled, p0i, p1i, 1, thickness=th)
                cv2.line(local_skel, p0i, p1i, 1, thickness=1)

        fy, fx = np.nonzero(local_filled)
        sy, sx = np.nonzero(local_skel)
        filled_idx = ((fy + y_min).astype(np.int64) * self.width) + (fx + x_min).astype(np.int64)
        skel_idx = ((sy + y_min).astype(np.int64) * self.width) + (sx + x_min).astype(np.int64)
        return filled_idx, skel_idx

    def _ensure_segment_pixel_cache(self, neuron_id: int, seg: Dict):
        if 'filled_idx' in seg and 'skeleton_idx' in seg:
            return
        f_idx, s_idx = self._rasterize_segment_indices(seg.get('points', []), seg.get('thickness', 1.0))
        seg['filled_idx'] = f_idx
        seg['skeleton_idx'] = s_idx

    def _refresh_global_masks_at_indices(self, dirty_idx):
        if dirty_idx is None:
            return
        dirty_idx = np.asarray(dirty_idx, dtype=np.int64)
        if dirty_idx.size == 0:
            return
        idxs = np.unique(dirty_idx)
        if idxs.size == 0:
            return

        # Compose by neuron draw order (ascending IDs overwrite previous), equivalent to max id present.
        filled_owner = np.zeros(idxs.shape, dtype=np.uint16)
        skel_owner = np.zeros(idxs.shape, dtype=np.uint16)
        for i, neuron in enumerate(self.neurons):
            nid = np.uint16(neuron.id_int)
            filled_owner[self._branch_filled_counts[i].reshape(-1)[idxs] > 0] = nid
            skel_owner[self._branch_skeleton_counts[i].reshape(-1)[idxs] > 0] = nid

        soma_owner = self.somas_mask.ravel()[idxs]
        soma_skel_owner = self.soma_skeleton_mask.ravel()[idxs]
        rows, cols = np.divmod(idxs, self.width)   # 2-D writes: safe even for non-contiguous masks
        self.filled_mask[rows, cols] = np.where(filled_owner > 0, filled_owner, soma_owner)
        self.skeleton_mask[rows, cols] = np.where(skel_owner > 0, skel_owner, soma_skel_owner)

    def _apply_segment_counts(self, neuron_id: int, seg: Dict, sign: int):
        """Add (+1) or remove (-1) one segment's pixels from the per-neuron count rasters.
        Returns the touched flat pixel indices. Indices within a segment are unique
        (np.nonzero output), so plain fancy-index += is exact (no np.add.at needed)."""
        self._ensure_segment_pixel_cache(neuron_id, seg)
        idx = int(neuron_id) - 1
        f_idx = seg.get('filled_idx', np.array([], dtype=np.int64))
        s_idx = seg.get('skeleton_idx', np.array([], dtype=np.int64))
        for counts, pix in ((self._branch_filled_counts[idx], f_idx),
                            (self._branch_skeleton_counts[idx], s_idx)):
            if pix.size:
                flat = counts.reshape(-1)
                if sign > 0:
                    flat[pix] += 1
                else:
                    flat[pix] = np.maximum(flat[pix] - 1, 0)
        return f_idx, s_idx

    def incremental_batch_update(self, neuron_id: int, removed_segs, added_segs):
        """Apply many segment removals/additions for one neuron with a single mask refresh."""
        if not self.incremental_branch_mode:
            return
        dirty = []
        for seg in removed_segs:
            dirty.extend(self._apply_segment_counts(neuron_id, seg, -1))
        for seg in added_segs:
            dirty.extend(self._apply_segment_counts(neuron_id, seg, +1))
        dirty = [d for d in dirty if d.size]
        if dirty:
            self._refresh_global_masks_at_indices(np.concatenate(dirty))

    def incremental_add_segment_record(self, neuron_id: int, seg: Dict):
        if not self.incremental_branch_mode:
            return
        f_idx, s_idx = self._apply_segment_counts(neuron_id, seg, +1)
        if f_idx.size or s_idx.size:
            self._refresh_global_masks_at_indices(np.concatenate([f_idx, s_idx]))

    def incremental_remove_segment_record(self, neuron_id: int, seg: Dict):
        if not self.incremental_branch_mode:
            return
        f_idx, s_idx = self._apply_segment_counts(neuron_id, seg, -1)
        if f_idx.size or s_idx.size:
            self._refresh_global_masks_at_indices(np.concatenate([f_idx, s_idx]))

    def rebuild_masks_from_branch_history(self, finalize_masks: bool=True,
                                         clear_interaction_annotations: bool=False):
        """
        Rebuild filled/skeleton masks from recorded branch segments, respecting each
        neuron's visible_total_length cutoff.
        """
        self.filled_mask.fill(0)
        self.skeleton_mask.fill(0)
        self.somas_mask.fill(0)
        replay_stats = {
            'segments_replayed': 0,
            'segments_clipped': 0,
            'segments_skipped': 0,
        }

        for neuron in self.neurons:
            # Redraw soma masks first (immutable)
            f, s = neuron.dendrite.soma.binary_masks((self.height, self.width), neuron.id_int)
            self.somas_mask[f > 0] = f[f > 0]
            self.filled_mask[f > 0] = f[f > 0]
            self.skeleton_mask[s > 0] = s[s > 0]

        for neuron in self.neurons:
            cutoff = float(getattr(neuron, 'visible_total_length', neuron.dendrite.total_length))
            for seg in neuron.branch_data:
                path_start = float(seg.get('path_start_dist', 0.0))
                path_end = float(seg.get('path_end_dist', path_start))
                if cutoff <= path_start + 1e-9:
                    replay_stats['segments_skipped'] += 1
                    continue

                points = seg.get('points', [])
                if not points or len(points) < 2:
                    replay_stats['segments_skipped'] += 1
                    continue

                if path_end <= path_start + 1e-9 or cutoff >= path_end - 1e-9:
                    draw_points = points
                else:
                    frac = (cutoff - path_start) / max(1e-9, (path_end - path_start))
                    draw_points = self._clip_polyline_to_fraction(points, frac)
                    if len(draw_points) < 2:
                        replay_stats['segments_skipped'] += 1
                        continue
                    replay_stats['segments_clipped'] += 1

                self._draw_history_polyline(neuron.id_int, draw_points, seg.get('thickness', 1.0))
                replay_stats['segments_replayed'] += 1

        if clear_interaction_annotations:
            for neuron in self.neurons:
                neuron.connections = []
                neuron.fusion_points = []
                neuron.soma_connections = []

        if finalize_masks:
            self._finalize_masks()
        self.last_rebuild_history_stats = replay_stats

    def apply_visible_lengths(self, total_lengths, finalize_masks: bool=True,
                              clear_interaction_annotations: bool=True):
        """
        Apply in-place prune/regrowth by changing visible length cutoffs and replaying
        recorded branch history into the masks.
        """
        if len(total_lengths) != len(self.neurons):
            raise ValueError("total_lengths length must match number of neurons")

        for neuron, length in zip(self.neurons, total_lengths):
            neuron.set_visible_total_length(float(length))

        self.rebuild_masks_from_branch_history(
            finalize_masks=finalize_masks,
            clear_interaction_annotations=clear_interaction_annotations
        )

    def _ids_present_near(self, x: float, y: float, radius: int) -> set:
        """Neuron ids whose arbor (dendrite count rasters, or soma) covers any pixel within
        `radius` of (x, y). Uses the per-neuron count layers in incremental mode, so
        overlapping arbors are all seen (the composed label mask only shows the max id)."""
        c, r = int(round(x)), int(round(y))
        r0, r1 = max(0, r - radius), min(self.height, r + radius + 1)
        c0, c1 = max(0, c - radius), min(self.width, c + radius + 1)
        if r0 >= r1 or c0 >= c1:
            return set()
        ids = set(int(v) for v in np.unique(self.somas_mask[r0:r1, c0:c1]) if v)
        if getattr(self, 'incremental_branch_mode', False):
            for i, neuron in enumerate(self.neurons):
                if np.any(self._branch_filled_counts[i][r0:r1, c0:c1] > 0):
                    ids.add(int(neuron.id_int))
        else:
            ids |= set(int(v) for v in np.unique(self.filled_mask[r0:r1, c0:c1]) if v)
        return ids

    def prune_stale_interactions(self, radius: int = 2) -> int:
        """
        Drop recorded interaction points (competition stops, fusions, soma contacts) whose
        geometry no longer exists after pruning/retraction: the neuron's own arbor must still
        reach the point, and for dendrite-dendrite interactions another neuron's arbor must
        still be there too. Returns the number of points removed.
        """
        removed = 0
        for neuron in self.neurons:
            nid = int(neuron.id_int)

            def alive(x, y, other=None, need_other=True):
                ids = self._ids_present_near(x, y, radius)
                if nid not in ids:
                    return False
                if not need_other:
                    return True
                if other is not None:
                    return int(other) in ids
                return len(ids - {nid}) > 0

            before = len(neuron.connections) + len(neuron.fusion_points) + len(neuron.soma_connections)
            neuron.connections = [p for p in neuron.connections if alive(p[0], p[1])]
            neuron.fusion_points = [p for p in neuron.fusion_points if alive(p[0], p[1], other=p[2])]
            # Soma contacts: somas never move, so only the dendrite end needs to survive.
            neuron.soma_connections = [p for p in neuron.soma_connections
                                       if alive(p[0], p[1], need_other=False)]
            removed += before - (len(neuron.connections) + len(neuron.fusion_points) + len(neuron.soma_connections))
        return removed

    def prune_and_forget_lengths(self, total_lengths, finalize_masks: bool=True,
                                 clear_interaction_annotations: bool=False):
        """
        Permanently prune recorded branch history beyond per-neuron cutoffs and
        reactivate frontier tips so subsequent regrowth can generate new paths.
        """
        if len(total_lengths) != len(self.neurons):
            raise ValueError("total_lengths length must match number of neurons")

        for neuron, length in zip(self.neurons, total_lengths):
            neuron.prune_and_forget(float(length))

        if self.incremental_branch_mode and not finalize_masks and not clear_interaction_annotations:
            return

        self.rebuild_masks_from_branch_history(
            finalize_masks=finalize_masks,
            clear_interaction_annotations=clear_interaction_annotations
        )

    def _apply_local_smoothing(self):
        """Apply local morphological smoothing around connection and fusion points."""
        if not np.any(self.filled_mask):
            return
            
        # Get all interaction points (connections and fusions)
        interaction_points = []
        
        # Collect connection points
        for neuron in self.neurons:
            for conn_x, conn_y, step in neuron.connections:
                interaction_points.append((int(conn_x), int(conn_y)))
                
        # Collect fusion points  
        for neuron in self.neurons:
            for fusion_x, fusion_y, other_id, step in neuron.fusion_points:
                interaction_points.append((int(fusion_x), int(fusion_y)))
        
        # Collect soma-dendrite connection points
        for neuron in self.neurons:
            for soma_conn_x, soma_conn_y, soma_id, step in neuron.soma_connections:
                interaction_points.append((int(soma_conn_x), int(soma_conn_y)))
        
        if not interaction_points:
            return
            
        # Create a mask for regions to smooth
        smooth_mask = np.zeros_like(self.filled_mask, dtype=bool)
        
        # Mark areas around interaction points for smoothing
        for x, y in interaction_points:
            x_min = max(0, x - self.smoothing_zone_radius)
            x_max = min(self.width, x + self.smoothing_zone_radius + 1)
            y_min = max(0, y - self.smoothing_zone_radius)
            y_max = min(self.height, y + self.smoothing_zone_radius + 1)
            
            # Create circular smoothing zone
            for yi in range(y_min, y_max):
                for xi in range(x_min, x_max):
                    if (xi - x)**2 + (yi - y)**2 <= self.smoothing_zone_radius**2:
                        smooth_mask[yi, xi] = True
        
        # Apply the selected smoothing method locally
        if self.local_smoothing_method == 'closing':
            self._apply_local_closing(smooth_mask)
        elif self.local_smoothing_method == 'gaussian':
            self._apply_local_gaussian(smooth_mask)
        elif self.local_smoothing_method == 'opening':
            self._apply_local_opening(smooth_mask)
        elif self.local_smoothing_method == 'distance':
            self._apply_local_distance_smoothing(smooth_mask)
        else:
            print(f"Warning: Unknown smoothing method '{self.local_smoothing_method}', using closing")
            self._apply_local_closing(smooth_mask)

    def _apply_local_closing(self, smooth_mask):
        """Apply morphological closing locally to smooth regions."""
        # Create structuring element
        selem = disk(self.smoothing_radius)
        
        # Apply closing only to the smoothing regions
        binary_mask = (self.filled_mask > 0)
        smoothed_binary = binary_closing(binary_mask, selem)
        
        # Only update regions marked for smoothing
        updated_mask = binary_mask.copy()
        updated_mask[smooth_mask] = smoothed_binary[smooth_mask]
        
        # Update filled_mask while preserving neuron IDs
        for neuron in self.neurons:
            neuron_mask = (self.filled_mask == neuron.id_int)
            if np.any(neuron_mask):
                # Apply the smoothing to this neuron's regions
                smoothed_neuron = neuron_mask.copy()
                smoothed_neuron[smooth_mask] = updated_mask[smooth_mask] & \
                    binary_closing(neuron_mask, selem)[smooth_mask]
                self.filled_mask[neuron_mask] = 0  # Clear old
                self.filled_mask[smoothed_neuron] = neuron.id_int  # Set new

    def _apply_local_opening(self, smooth_mask):
        """Apply morphological opening locally to smooth regions."""
        # Create structuring element
        selem = disk(self.smoothing_radius)
        
        # Apply opening only to the smoothing regions
        binary_mask = (self.filled_mask > 0)
        smoothed_binary = binary_opening(binary_mask, selem)
        
        # Only update regions marked for smoothing
        updated_mask = binary_mask.copy()
        updated_mask[smooth_mask] = smoothed_binary[smooth_mask]
        
        # Update filled_mask while preserving neuron IDs
        for neuron in self.neurons:
            neuron_mask = (self.filled_mask == neuron.id_int)
            if np.any(neuron_mask):
                # Apply the smoothing to this neuron's regions
                smoothed_neuron = neuron_mask.copy()
                smoothed_neuron[smooth_mask] = updated_mask[smooth_mask] & \
                    binary_opening(neuron_mask, selem)[smooth_mask]
                self.filled_mask[neuron_mask] = 0  # Clear old
                self.filled_mask[smoothed_neuron] = neuron.id_int  # Set new

    def _apply_local_gaussian(self, smooth_mask):
        """Apply Gaussian smoothing with thresholding locally."""
        # For each neuron, apply Gaussian smoothing to its regions
        for neuron in self.neurons:
            neuron_mask = (self.filled_mask == neuron.id_int).astype(float)
            if np.any(neuron_mask):
                # Apply Gaussian smoothing
                smoothed = gaussian(neuron_mask, sigma=self.gaussian_sigma)
                
                # Threshold back to binary (0.5 threshold)
                smoothed_binary = smoothed > 0.5
                
                # Only update regions marked for smoothing
                updated_neuron_mask = (self.filled_mask == neuron.id_int)
                updated_neuron_mask[smooth_mask] = smoothed_binary[smooth_mask]
                
                # Update the mask
                self.filled_mask[self.filled_mask == neuron.id_int] = 0  # Clear old
                self.filled_mask[updated_neuron_mask] = neuron.id_int  # Set new

    def _apply_local_distance_smoothing(self, smooth_mask):
        """Apply distance transform-based smoothing for very smooth transitions."""
        # For each neuron, apply distance transform smoothing
        for neuron in self.neurons:
            neuron_mask = (self.filled_mask == neuron.id_int)
            if np.any(neuron_mask):
                # Compute distance transform from the boundary
                dist_transform = distance_transform_edt(neuron_mask)
                
                # Create smooth transition using distance transform
                # Smoothing strength decreases with distance from interaction points
                smooth_strength = np.zeros_like(dist_transform)
                for y in range(self.height):
                    for x in range(self.width):
                        if smooth_mask[y, x]:
                            # Distance-based smoothing: closer to interaction points = more smoothing
                            smooth_strength[y, x] = 1.0
                
                # Apply smoothing by modifying the distance field
                smoothed_dist = gaussian(dist_transform * smooth_strength, sigma=self.gaussian_sigma)
                
                # Threshold back to binary using original boundary distance
                smoothed_binary = smoothed_dist > (self.smoothing_radius * 0.5)
                
                # Only update regions marked for smoothing  
                updated_neuron_mask = neuron_mask.copy()
                updated_neuron_mask[smooth_mask] = smoothed_binary[smooth_mask]
                
                # Update the mask
                self.filled_mask[neuron_mask] = 0  # Clear old
                self.filled_mask[updated_neuron_mask] = neuron.id_int  # Set new

    def get_all_connections(self):
        """Return all connection points across all neurons."""
        all_connections = []
        for neuron in self.neurons:
            for conn_point in neuron.connections:
                all_connections.append({
                    'neuron_id': neuron.nid_str,
                    'neuron_int_id': neuron.id_int,
                    'connection_point': conn_point
                })
        return all_connections

    def analyze_growth_statistics(self):
        """Return statistics about growth including connections, pass-through, and fusion behavior."""
        total_connections = sum(len(n.connections) for n in self.neurons)
        total_fusions = sum(len(n.fusion_points) for n in self.neurons)
        total_soma_connections = sum(len(n.soma_connections) for n in self.neurons)
        total_active_tips = sum(len(n.active_tips) for n in self.neurons)
        growing_neurons = sum(1 for n in self.neurons if n.is_growing)
        
        return {
            'total_connections': total_connections,
            'connections_per_neuron': [len(n.connections) for n in self.neurons],
            'total_fusions': total_fusions,
            'fusions_per_neuron': [len(n.fusion_points) for n in self.neurons],
            'total_soma_connections': total_soma_connections,
            'soma_connections_per_neuron': [len(n.soma_connections) for n in self.neurons],
            'total_active_tips': total_active_tips,
            'growing_neurons': growing_neurons,
            'average_pass_through_prob': sum(n.pass_prob for n in self.neurons) / len(self.neurons) if self.neurons else 0,
            'average_fusion_prob': sum(n.fusion_prob for n in self.neurons) / len(self.neurons) if self.neurons else 0
        }

    def render_frame(self, mask_type='filled', colormap='viridis', show_soma_labels=False, 
                    show_connections=False, show_fusions=False, show_soma_connections=False, figsize=(10, 10), dpi=100,
                    historical_connections=None, historical_fusions=None, historical_soma_connections=None,
                    background_color='black'):
        """
        Render a single frame of the network state.
        
        Args:
            mask_type: 'filled', 'skeleton', or 'somas'
            colormap: matplotlib colormap name
            show_soma_labels: whether to show neuron ID labels
            show_connections: whether to highlight connection points
            show_fusions: whether to highlight fusion points
            show_soma_connections: whether to highlight soma-dendrite connection points
            figsize: figure size in inches
            dpi: dots per inch for rendering
            background_color: 'black' or 'white' background color
            
        Returns:
            numpy array of the rendered image (RGB)
        """
        plt.ioff()  # Turn off interactive mode
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        
        # Select the appropriate mask
        if mask_type == 'filled':
            mask = self.filled_mask
        elif mask_type == 'skeleton':
            mask = self.skeleton_mask
        elif mask_type == 'somas':
            mask = self.somas_mask
        else:
            raise ValueError("mask_type must be 'filled', 'skeleton', or 'somas'")
        
        # Create colored visualization
        if mask.max() > 0:
            normalized_mask = mask.astype(float) / mask.max()
            colored_mask = plt.cm.get_cmap(colormap)(normalized_mask)
            # Set background color based on parameter
            if background_color.lower() == 'black':
                colored_mask[mask == 0] = [0, 0, 0, 1]  # Black background
            else:
                colored_mask[mask == 0] = [1, 1, 1, 1]  # White background
        else:
            # Set background when no data
            if background_color.lower() == 'black':
                colored_mask = np.zeros((*mask.shape, 4))  # Black background when no data
            else:
                colored_mask = np.ones((*mask.shape, 4))   # White background when no data
        
        ax.imshow(colored_mask, origin='upper')
        
        # Add soma labels if requested
        if show_soma_labels:
            for neuron in self.neurons:
                x, y = neuron.position[:2]
                ax.text(x, y, str(neuron.id_int), ha='center', va='center',
                       color='black', fontsize=8, fontweight='bold',
                       bbox=dict(boxstyle='circle,pad=0.2', facecolor='white', alpha=0.8))
        
        # Add connection points if requested
        if show_connections:
            if historical_connections is not None:
                # Use historical data (for growth history visualization)
                for conn_x, conn_y, step in historical_connections:
                    ax.plot(conn_x, conn_y, 'ro', markersize=3, alpha=0.7)
            else:
                # Use current neuron data (for final state visualization)
                for neuron in self.neurons:
                    for conn_x, conn_y, step in neuron.connections:
                        ax.plot(conn_x, conn_y, 'ro', markersize=3, alpha=0.7)
        
        # Add fusion points if requested
        if show_fusions:
            if historical_fusions is not None:
                # Use historical data (for growth history visualization)
                for fusion_x, fusion_y, other_id, step in historical_fusions:
                    ax.plot(fusion_x, fusion_y, 'go', markersize=4, alpha=0.8)  # Green circles for fusions
            else:
                # Use current neuron data (for final state visualization)
                for neuron in self.neurons:
                    for fusion_x, fusion_y, other_id, step in neuron.fusion_points:
                        ax.plot(fusion_x, fusion_y, 'go', markersize=4, alpha=0.8)  # Green circles for fusions
        
        # Add soma-dendrite connection points if requested
        if show_soma_connections:
            if historical_soma_connections is not None:
                # Use historical data (for growth history visualization)
                for soma_conn_x, soma_conn_y, soma_id, step in historical_soma_connections:
                    ax.plot(soma_conn_x, soma_conn_y, 'bo', markersize=3, alpha=0.7)
            else:
                # Use current neuron data (for final state visualization)
                for neuron in self.neurons:
                    for soma_conn_x, soma_conn_y, soma_id, step in neuron.soma_connections:
                        ax.plot(soma_conn_x, soma_conn_y, 'bo', markersize=3, alpha=0.7)
        
        ax.set_xlim(0, self.width)
        ax.set_ylim(self.height, 0)  # Flip y-axis to match image coordinates
        ax.set_aspect('equal')
        ax.axis('off')
        
        # Convert to numpy array
        fig.canvas.draw()
        buf = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
        buf = buf.reshape(fig.canvas.get_width_height()[::-1] + (3,))
        
        plt.close(fig)
        return buf
    
    def render_growth_history(self, output_dir='network_growth_frames', mask_type='filled', 
                             colormap='viridis', show_soma_labels=False, show_connections=False,
                             show_fusions=False, show_soma_connections=False, figsize=(10, 10), dpi=100,
                             background_color='black'):
        """
        Render all frames from the growth history.
        
        Args:
            output_dir: directory to save frame images
            mask_type: 'filled', 'skeleton', or 'somas'  
            colormap: matplotlib colormap name
            show_soma_labels: whether to show neuron ID labels
            show_connections: whether to highlight connection points
            show_fusions: whether to highlight fusion points
            show_soma_connections: whether to highlight soma-dendrite connection points
            figsize: figure size in inches
            dpi: dots per inch for rendering
            
        Returns:
            list of frame file paths
        """
        if not self.growth_history_masks:
            print("No growth history recorded. Set record_history_every_n_steps > 0 in sim_props.")
            return []
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        frame_paths = []
        
        # Save original state
        original_filled = self.filled_mask.copy()
        original_skeleton = self.skeleton_mask.copy()
        original_somas = self.somas_mask.copy()
        
        try:
            for i, history_entry in enumerate(self.growth_history_masks):
                # Temporarily replace masks with historical data
                self.filled_mask = history_entry['filled_mask']
                self.skeleton_mask = history_entry['skeleton_mask'] 
                self.somas_mask = history_entry['somas_mask']
                
                # Render frame
                frame_img = self.render_frame(
                    mask_type=mask_type,
                    colormap=colormap,
                    show_soma_labels=show_soma_labels,
                    show_connections=show_connections,
                    show_fusions=show_fusions,
                    show_soma_connections=show_soma_connections,
                    figsize=figsize,
                    dpi=dpi,
                    historical_connections=history_entry.get('connections', []),
                    historical_fusions=history_entry.get('fusions', []),
                    historical_soma_connections=history_entry.get('soma_connections', []),
                    background_color=background_color
                )
                
                # Save frame
                frame_path = os.path.join(output_dir, f'frame_{i:04d}_step_{history_entry["step"]:04d}.png')
                Image.fromarray(frame_img).save(frame_path)
                frame_paths.append(frame_path)
                
        finally:
            # Restore original state
            self.filled_mask = original_filled
            self.skeleton_mask = original_skeleton
            self.somas_mask = original_somas
        
        print(f"Rendered {len(frame_paths)} frames in {output_dir}")
        return frame_paths
    
    def create_growth_gif(self, gif_path='network_growth.gif', mask_type='filled',
                         colormap='viridis', show_soma_labels=False, show_connections=False,
                         show_fusions=False, show_soma_connections=False, figsize=(10, 10), dpi=100, duration=200, loop=0,
                         background_color='black'):
        """
        Create an animated GIF of the network growth.
        
        Args:
            gif_path: path to save the GIF file
            mask_type: 'filled', 'skeleton', or 'somas'
            colormap: matplotlib colormap name
            show_soma_labels: whether to show neuron ID labels
            show_connections: whether to highlight connection points
            show_fusions: whether to highlight fusion points
            show_soma_connections: whether to highlight soma-dendrite connection points
            figsize: figure size in inches
            dpi: dots per inch for rendering
            duration: milliseconds between frames
            loop: number of loops (0 = infinite)
            
        Returns:
            path to the created GIF file
        """
        if not self.growth_history_masks:
            print("No growth history recorded. Set record_history_every_n_steps > 0 in sim_props.")
            return None
        
        frames = []
        
        # Save original state
        original_filled = self.filled_mask.copy()
        original_skeleton = self.skeleton_mask.copy()
        original_somas = self.somas_mask.copy()
        
        try:
            print(f"Creating GIF with {len(self.growth_history_masks)} frames...")
            
            for i, history_entry in enumerate(self.growth_history_masks):
                # Temporarily replace masks with historical data
                self.filled_mask = history_entry['filled_mask']
                self.skeleton_mask = history_entry['skeleton_mask']
                self.somas_mask = history_entry['somas_mask']
                
                # Render frame
                frame_img = self.render_frame(
                    mask_type=mask_type,
                    colormap=colormap,
                    show_soma_labels=show_soma_labels,
                    show_connections=show_connections,
                    show_fusions=show_fusions,
                    show_soma_connections=show_soma_connections,
                    figsize=figsize,
                    dpi=dpi,
                    historical_connections=history_entry.get('connections', []),
                    historical_fusions=history_entry.get('fusions', []),
                    historical_soma_connections=history_entry.get('soma_connections', []),
                    background_color=background_color
                )
                
                frames.append(Image.fromarray(frame_img))
                
                if (i + 1) % 10 == 0:
                    print(f"Processed {i + 1}/{len(self.growth_history_masks)} frames")
                    
        finally:
            # Restore original state
            self.filled_mask = original_filled
            self.skeleton_mask = original_skeleton
            self.somas_mask = original_somas
        
        # Save GIF
        if frames:
            frames[0].save(
                gif_path,
                save_all=True,
                append_images=frames[1:],
                duration=duration,
                loop=loop,
                optimize=True
            )
            print(f"GIF saved to: {gif_path}")
            return gif_path
        else:
            print("No frames generated.")
            return None
    
    def create_comparison_gif(self, gif_path='network_comparison.gif', figsize=(12, 5), 
                             dpi=100, duration=200, loop=0, background_color='black'):
        """
        Create a GIF showing filled and skeleton views side by side.
        
        Args:
            gif_path: path to save the GIF file
            figsize: figure size in inches
            dpi: dots per inch for rendering
            duration: milliseconds between frames
            loop: number of loops (0 = infinite)
            background_color: 'black' or 'white' background color
            
        Returns:
            path to the created GIF file
        """
        if not self.growth_history_masks:
            print("No growth history recorded. Set record_history_every_n_steps > 0 in sim_props.")
            return None
        
        frames = []
        
        # Save original state
        original_filled = self.filled_mask.copy()
        original_skeleton = self.skeleton_mask.copy()
        original_somas = self.somas_mask.copy()
        
        try:
            print(f"Creating comparison GIF with {len(self.growth_history_masks)} frames...")
            
            for i, history_entry in enumerate(self.growth_history_masks):
                plt.ioff()
                fig, axes = plt.subplots(1, 2, figsize=figsize, dpi=dpi)
                
                # Temporarily replace masks with historical data
                filled = history_entry['filled_mask']
                skeleton = history_entry['skeleton_mask']
                
                # Render each view with summer colormap
                for ax, mask, title in zip(
                    axes, 
                    [filled, skeleton],
                    ['Filled Dendrites', 'Skeleton']
                ):
                    if mask.max() > 0:
                        normalized_mask = mask.astype(float) / mask.max()
                        colored_mask = plt.cm.get_cmap('summer')(normalized_mask)
                        # Set background color based on parameter
                        if background_color.lower() == 'black':
                            colored_mask[mask == 0] = [0, 0, 0, 1]  # Black background
                        else:
                            colored_mask[mask == 0] = [1, 1, 1, 1]  # White background
                    else:
                        # Set background when no data
                        if background_color.lower() == 'black':
                            colored_mask = np.zeros((*mask.shape, 4))  # Black background when no data
                        else:
                            colored_mask = np.ones((*mask.shape, 4))   # White background when no data
                    
                    ax.imshow(colored_mask, origin='upper')
                    ax.set_title(f'{title} (Step {history_entry["step"]})', fontsize=10)
                    ax.set_xlim(0, self.width)
                    ax.set_ylim(self.height, 0)
                    ax.set_aspect('equal')
                    ax.axis('off')
                
                plt.tight_layout()
                
                # Convert to numpy array
                fig.canvas.draw()
                buf = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
                buf = buf.reshape(fig.canvas.get_width_height()[::-1] + (3,))
                frames.append(Image.fromarray(buf))
                
                plt.close(fig)
                
                if (i + 1) % 10 == 0:
                    print(f"Processed {i + 1}/{len(self.growth_history_masks)} frames")
                    
        finally:
            # Restore original state
            self.filled_mask = original_filled
            self.skeleton_mask = original_skeleton
            self.somas_mask = original_somas
        
        # Save GIF
        if frames:
            frames[0].save(
                gif_path,
                save_all=True,
                append_images=frames[1:],
                duration=duration,
                loop=loop,
                optimize=True
            )
            print(f"Comparison GIF saved to: {gif_path}")
            return gif_path
        else:
            print("No frames generated.")
            return None

################################################################################
#                     EXPOSED HIGH‑LEVEL HELPER FUNCTION                        #
################################################################################

def generate_network(network_params:dict, neuron_params:dict, sim_props:dict=None, 
                     simage_props:dict=None, network_id:str=None, seed_coordinates=None):
    
    # Extract basic network parameters
    width = network_params['width']
    height = network_params['height']
    num_neurons = network_params['num_neurons']
    edge_margin = network_params.get('edge_margin', 50)
    deterministic = network_params.get('deterministic', False)
    
    # Create network
    network = Network(
        width=width, 
        height=height, 
        num_neurons=num_neurons,
        neuron_params=neuron_params,
        edge_margin=edge_margin,
        sim_props=sim_props,
        deterministic=deterministic,
        seed_coordinates=seed_coordinates,
        network_id=network_id,
        simage_props=simage_props,
        network_params=network_params
    )
    network.grow_network()
    return network
