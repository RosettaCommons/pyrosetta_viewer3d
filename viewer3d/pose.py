import attr
import logging

from pyrosetta import Pose
from typing import List, NoReturn, Optional

from viewer3d.validators import requires_show


_logger: logging.Logger = logging.getLogger("viewer3d.pose")


@attr.s(kw_only=False, slots=False)
class PoseBase:
    @requires_show
    def add_pose(
        self, pose: Pose, index: Optional[int] = None, update_viewer: bool = True
    ) -> None:
        if index is None:
            index = self.get_decoy_widget_index()
        self.poses[index].append(pose)
        self.pdbstrings[index].append(None)
        if update_viewer:
            model = self.get_last_model(index)
            self.update_viewer(
                index=index, model=model, add_objects=True, remove_objects=False
            )

    @requires_show
    def add_pdbstring(
        self, pdbstring: str, index: Optional[int] = None, update_viewer: bool = True
    ) -> None:
        if index is None:
            index = self.get_decoy_widget_index()
        self.poses[index].append(None)
        self.pdbstrings[index].append(pdbstring)
        if update_viewer:
            model = self.get_last_model(index)
            self.update_viewer(
                index=index, model=model, add_objects=True, remove_objects=False
            )

    def remove_pose(
        self,
        index: Optional[int] = None,
        model: Optional[int] = None,
        update_viewer: bool = True,
    ) -> Optional[NoReturn]:
        self.remove_pdbstring(index=index, model=model, update_viewer=update_viewer)

    @requires_show
    def remove_pdbstring(
        self,
        index: Optional[int] = None,
        model: Optional[int] = None,
        update_viewer: bool = True,
    ) -> Optional[NoReturn]:
        if index is None:
            index = self.get_decoy_widget_index()
        if model is None or model not in set(range(len(self.pdbstrings[index]))):
            model = self.get_last_model(index)
        if len(self.pdbstrings[index]) > 0:
            self.poses[index].pop(model)
            self.pdbstrings[index].pop(model)
        else:
            raise IndexError(
                f"The 'poses' and 'pdbstrings' attributes are empty at index `{index}`."
            )
        if update_viewer:
            self.update_viewer(
                index=index, model=model, add_objects=False, remove_objects=True
            )

    @requires_show
    def update_pose(
        self,
        pose: Pose,
        index: Optional[int] = None,
        model: Optional[int] = None,
        update_viewer: bool = True,
    ) -> Optional[NoReturn]:
        if index is None:
            index = self.get_decoy_widget_index()
        if model is None or model not in set(range(len(self.poses[index]))):
            model = 0
        if index in self.poses.keys():
            self.poses[index][model] = pose
            self.pdbstrings[index][model] = None
        else:
            raise IndexError(
                f"The 'poses' and 'pdbstrings' attributes do not have index `{index}`."
            )
        if update_viewer:
            self.update_viewer(
                index=index, model=model, add_objects=True, remove_objects=True
            )

    @requires_show
    def update_pdbstring(
        self,
        pdbstring: str,
        index: Optional[int] = None,
        model: Optional[int] = None,
        update_viewer: bool = True,
    ) -> Optional[NoReturn]:
        if index is None:
            index = self.get_decoy_widget_index()
        if model is None or model not in set(range(len(self.pdbstrings[index]))):
            model = 0
        if index in self.pdbstrings.keys():
            self.poses[index][model] = None
            self.pdbstrings[index][model] = pdbstring
        else:
            raise IndexError(
                f"The 'poses' and 'pdbstrings' attributes do not have index `{index}`."
            )
        if update_viewer:
            self.update_viewer(
                index=index, model=model, add_objects=True, remove_objects=True
            )

    @requires_show
    def update_poses(
        self, poses: List[Pose], index: Optional[int] = None, update_viewer: bool = True
    ) -> Optional[NoReturn]:
        if index is None:
            index = self.get_decoy_widget_index()
        assert isinstance(poses, list), "The 'poses' argument must be of type `list`."
        for pose in poses:
            assert isinstance(
                pose, Pose
            ), f"Object must be of type `Pose`. Received: {type(pose)}"
        self.poses[index] = poses
        self.pdbstrings[index] = [None] * len(poses)
        if update_viewer:
            self.update_viewer(
                index=index, model=None, add_objects=True, remove_objects=True
            )

    @requires_show
    def update_pdbstrings(
        self,
        pdbstrings: List[str],
        index: Optional[int] = None,
        update_viewer: bool = True,
    ) -> Optional[NoReturn]:
        if index is None:
            index = self.get_decoy_widget_index()
        assert isinstance(
            pdbstrings, list
        ), "The 'pdbstrings' argument must be of type `list`."
        for pdbstring in pdbstrings:
            assert isinstance(
                pdbstring, str
            ), f"Object must be of type `str`. Received: {type(pdbstring)}"
        self.poses[index] = [None] * len(pdbstrings)
        self.pdbstrings[index] = pdbstrings
        if update_viewer:
            self.update_viewer(
                index=index, model=None, add_objects=True, remove_objects=True
            )

    @requires_show
    def overlay(self, update_viewer: bool = True) -> None:
        """Re-index poses and pdbstrings into a single list."""
        for index in sorted(self.poses.keys()):
            if index != 0:
                self.poses[0].extend(self.poses[index])
                self.pdbstrings[0].extend(self.pdbstrings[index])
                self.poses.pop(index)
                self.pdbstrings.pop(index)
        if update_viewer:
            self.update_viewer(index=0)
