from fastapi import HTTPException
from services.informatic_service import InformaticsService, InformaticsCreate, Informatics

class InformaticsController:
    def __init__(self):
        self.informatic_service = InformaticsService()

    def create_informatic(self, informatic: InformaticsCreate):
        return self.informatic_service.create_informatic(informatic)

    def get_informatics(self):
        return self.informatic_service.get_informatics()

    def get_informatic(self, informatic_id: int):
        informatic = self.informatic_service.get_informatic(informatic_id)
        if informatic is None:
            raise HTTPException(status_code=404, detail=f"Informatic resource #{informatic_id} not found")
        return informatic

    def update_informatic(self, informatic_id: int, updated_informatic: InformaticsCreate):
        informatic = self.informatic_service.update_informatic(informatic_id, updated_informatic)
        if informatic is None:
            raise HTTPException(status_code=404, detail=f"Informatic resource #{informatic_id} not found")
        return informatic

    def delete_informatic(self, informatic_id: int):
        if self.informatic_service.delete_informatic(informatic_id):
            return {"message": f"Informatic resource #{informatic_id} deleted successfully"}
        raise HTTPException(status_code=404, detail=f"Informatic resource #{informatic_id} not found")