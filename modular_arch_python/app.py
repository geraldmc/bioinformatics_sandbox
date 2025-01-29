from fastapi import FastAPI
from routers import informatic_router

app = FastAPI()

app.include_router(informatic_router.router)

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("app:app", port=8000, host="127.0.0.1", reload=True)