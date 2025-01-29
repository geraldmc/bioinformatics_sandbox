from fastapi import APIRouter
from controllers.informatic_controller import InformaticsController, InformaticsCreate, Informatics

router = APIRouter()
informatic_controller = InformaticsController()

@router.post("/informatics", response_model=Informatics)
def create_informatic(informatic: InformaticsCreate):
  return informatic_controller.create_informatic(informatic)

@router.get("/informatics", response_model=list[Informatics])
def get_informatics():
  return informatic_controller.get_informatics()

@router.get("/informatics/{informatic_id}", response_model=Informatics)
def get_informatic(informatic_id: int):
  return informatic_controller.get_informatic(informatic_id)

@router.put("/informatics/{informatic_id}", response_model=Informatics)
def update_informatic(informatic_id: int, updated_informatic: InformaticsCreate):
  return informatic_controller.update_informatic(informatic_id, updated_informatic)

@router.delete("/informatics/{informatic_id}")
def delete_informatic(informatic_id: int):
  return informatic_controller.delete_informatic(informatic_id)